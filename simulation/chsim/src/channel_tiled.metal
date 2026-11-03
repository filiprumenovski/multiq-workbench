// Tiled / batched execution of the channel kernel. Compiled AFTER channel.metal in the same translation unit
// (host concatenates "#define ..." + channel.metal + channel_tiled.metal), so ChParams, ChJob, ChIon, rand01, gauss3,
// iso_dir are shared with the reference kernel `channel`. Every kernel here performs the same float arithmetic per step as
// `channel` (verified bitwise by `chsim check-tiled`); only WHERE the field values come from and HOW ions are scheduled differ.
//
// Field layout ("bricks"): the grid is cut into BR=7-cell bricks; brick k holds nodes 7k..7k+7 (8 per axis, shared faces
// duplicated), each node as three consecutive float4 {RF (w = conductor mask), combined DC, AC}. Every interpolation cell lies
// inside exactly one brick, so a lookup touches one 24 KB brick (1 page, a few cache lines) instead of 3 arrays x 4 rows spread
// over MBs, and a threadgroup can stage a whole block with one contiguous copy.
//
// Execution:
//   * Live ions are kept in a dense compacted list (tc_filter / tc_count), so finished ions never occupy lanes.
//   * tc_tail: persistent workers pull ions from that list and advance each by `chunk` steps with brick lookups (direct path).
//     Register-cached interpolation corners, next-step field/mask loads and RF/AC drive issued before the collision code.
//   * Optional tiled rounds (CHSIM_ROUNDS=<min live>): tc_clear -> tc_count (list compaction + block rank) -> tc_scan -> tc_scatter
//     -> tc_push, one threadgroup per (block, sub-batch of TILE_TG ions); the block's 512 nodes are staged in threadgroup memory,
//     each thread advances its ion up to `chunk` steps and parks (state untouched) when its cell leaves the block.

#ifndef TILE_TG
#define TILE_TG 256
#endif
#ifndef TAIL_TG
#define TAIL_TG 128
#endif
#ifndef TAIL_WMAX
#define TAIL_WMAX 16384
#endif
#ifndef TAIL_LSPAN
#define TAIL_LSPAN 64     // tc_size_args: lanes per simdgroup L = smallest power of two with count <= LSPAN * L (capped at 32)
#endif
#ifndef CH_NEAREST
#define CH_NEAREST 0
#endif

#define BR 7
#define BN 8
#define BNODES 512

struct TParams {
    int nbx, nby, nbz, nblocks;
    int chunk, n_ions, pad0, pad1;
};

// float4 index of node n's RF entry (DC at +1, AC at +2) inside the brick layout.
static inline int brick_node(int3 n, int nbx, int nby) {
    int3 bi = n / BR; int3 li = n - bi * BR;
    return ((bi.z * nby + bi.y) * nbx + bi.x) * (BNODES * 3) + 3 * (li.x + BN * (li.y + BN * li.z));
}

#define NBX ((P.nx - 1 + BR - 1) / BR)
#define NBY ((P.ny - 1 + BR - 1) / BR)

// Build the brick layout from the three separate field arrays (once per DC configuration).
kernel void tc_brick(device const float4* Frf [[buffer(1)]], device const float4* Fdc [[buffer(2)]],
                     device const float4* Fac [[buffer(3)]], constant ChParams& P [[buffer(4)]],
                     device float4* out [[buffer(22)]], uint e [[thread_position_in_grid]])
{
    const uint brick = e / BNODES, ln = e % BNODES;
    const int lx = ln % BN, ly = (ln / BN) % BN, lz = ln / (BN * BN);
    const int nbx = NBX, nby = NBY;
    const int bx = brick % nbx, by = (brick / nbx) % nby, bz = brick / (nbx * nby);
    const int gx = min(bx * BR + lx, P.nx - 1), gy = min(by * BR + ly, P.ny - 1), gz = min(bz * BR + lz, P.nz - 1);
    const int gi = gx + P.nx * (gy + P.ny * gz);
    float4 acv = Fac[gi];
#ifndef AC_FULL
    if (P.nearest == 0) {   // AC slot of the brick = average over the 8 corners of the cell whose lower corner is this node (see CellCache)
        float3 a = float3(0.0f);
        for (int q = 0; q < 8; ++q) {
            int3 g = min(int3(gx + (q & 1), gy + ((q >> 1) & 1), gz + (q >> 2)), int3(P.nx - 1, P.ny - 1, P.nz - 1));
            a += Fac[g.x + P.nx * (g.y + P.ny * g.z)].xyz;
        }
        acv = float4(a * 0.125f, acv.w);
    }
#endif
    out[3 * e] = Frf[gi]; out[3 * e + 1] = Fdc[gi]; out[3 * e + 2] = acv;
}

static inline float3 tc_gc(int3 p, constant ChParams& P) {
    const float inv_scale = 1.0f / P.pos_scale;
    float3 x = float3(p) * inv_scale;
    return clamp(x * P.inv_h, float3(0.0f), float3(P.nx - 1.001f, P.ny - 1.001f, P.nz - 1.001f));
}

// Trilinear interpolation at node index b with node stride S (3 in the global brick layout where k picks the basis, 1 in a
// threadgroup array holding one basis) -- same corner order and mix() chain as `interp` in channel.metal.
template <int S, typename Ptr>
static inline float3 samp_tri(Ptr F, int b, int k, float3 t) {
    float3 c000 = F[b + k].xyz,                      c100 = F[b + S + k].xyz;
    float3 c010 = F[b + S * BN + k].xyz,             c110 = F[b + S * BN + S + k].xyz;
    float3 c001 = F[b + S * BN * BN + k].xyz,        c101 = F[b + S * BN * BN + S + k].xyz;
    float3 c011 = F[b + S * (BN + BN * BN) + k].xyz, c111 = F[b + S * (BN + BN * BN) + S + k].xyz;
    float3 c00 = mix(c000, c100, t.x), c10 = mix(c010, c110, t.x);
    float3 c01 = mix(c001, c101, t.x), c11 = mix(c011, c111, t.x);
    return mix(mix(c00, c10, t.y), mix(c01, c11, t.y), t.z);
}

// Register-cached corner data: an ion changes interpolation cell rarely (S14 ions move ~0.04 cell/step), so the 8 corners x 3 bases
// are fetched only when the cell index changes and the mix() chain runs from registers. Values and arithmetic are unchanged.
#ifdef AC_FULL
struct CellCache { float3 c[3][8]; int3 cell; };           // reference: RF, DC and AC corners all register-cached (72 floats), bitwise = direct kernel
#else
// Default: RF and DC corners cached as named scalars (arrays here were demoted to stack memory by the compiler), the AC field cached as its
// cell average (3 floats instead of 24; the brick's AC slot already holds that average, see tc_brick). The AC dipole field is uniform to
// ~3e-5 over r0 in the bulk and varies <= 10% per cell only in the entrance fringe where its force is ~500x below the RF force.
// Frees ~21 registers (JUMP.md section 1d). `CHSIM_DEFS=AC_FULL` restores the bitwise-reference arithmetic.
struct CellCache { float3 r0, r1, r2, r3, r4, r5, r6, r7, d0, d1, d2, d3, d4, d5, d6, d7; float3 ac; int3 cell; };
#endif
template <typename Ptr>
static inline void cell_load(thread CellCache& cc, Ptr F, int b, int S, int stride_basis, int3 i) {
    cc.cell = i;
    const int o1 = S, o2 = S * BN, o3 = S * BN + S, o4 = S * BN * BN, o5 = S * BN * BN + S, o6 = S * (BN + BN * BN), o7 = S * (BN + BN * BN) + S;
#ifdef AC_FULL
    const int off[8] = { 0, o1, o2, o3, o4, o5, o6, o7 };
#endif
#ifdef AC_FULL
    for (int k = 0; k < 3; ++k) for (int q = 0; q < 8; ++q) cc.c[k][q] = F[b + off[q] + k * stride_basis].xyz;
#else
    cc.r0 = F[b].xyz; cc.r1 = F[b + o1].xyz; cc.r2 = F[b + o2].xyz; cc.r3 = F[b + o3].xyz;
    cc.r4 = F[b + o4].xyz; cc.r5 = F[b + o5].xyz; cc.r6 = F[b + o6].xyz; cc.r7 = F[b + o7].xyz;
    cc.d0 = F[b + stride_basis].xyz; cc.d1 = F[b + o1 + stride_basis].xyz; cc.d2 = F[b + o2 + stride_basis].xyz; cc.d3 = F[b + o3 + stride_basis].xyz;
    cc.d4 = F[b + o4 + stride_basis].xyz; cc.d5 = F[b + o5 + stride_basis].xyz; cc.d6 = F[b + o6 + stride_basis].xyz; cc.d7 = F[b + o7 + stride_basis].xyz;
    cc.ac = F[b + 2 * stride_basis].xyz;               // brick AC slot = cell average (built by tc_brick)
#endif
}
static inline float3 mix8(float3 a0, float3 a1, float3 a2, float3 a3, float3 a4, float3 a5, float3 a6, float3 a7, float3 t) {
    float3 c00 = mix(a0, a1, t.x), c10 = mix(a2, a3, t.x);
    float3 c01 = mix(a4, a5, t.x), c11 = mix(a6, a7, t.x);
    return mix(mix(c00, c10, t.y), mix(c01, c11, t.y), t.z);
}
static inline float3 cell_mix(thread const float3* c, float3 t) {
    float3 c00 = mix(c[0], c[1], t.x), c10 = mix(c[2], c[3], t.x);
    float3 c01 = mix(c[4], c[5], t.x), c11 = mix(c[6], c[7], t.x);
    return mix(mix(c00, c10, t.y), mix(c01, c11, t.y), t.z);
}


// Same arithmetic as gauss3/iso_dir in channel.metal, with sin/cos of the common angle evaluated as one sincos() (bitwise equal to
// separate sin()/cos(); verified on 1e6 angles and by check-tiled).
static inline float3 gauss3_sc(thread uint4 &st) {
    float u1 = max(rand01(st), 1e-7f), u2 = rand01(st), u3 = max(rand01(st), 1e-7f), u4 = rand01(st);
    float r1 = sqrt(-2.0f * log(u1)), r2 = sqrt(-2.0f * log(u3));
    float c1, s1 = sincos(6.2831853f * u2, c1);
    float c2 = cos(6.2831853f * u4);
    return float3(r1 * c1, r1 * s1, r2 * c2);
}
static inline float3 iso_dir_sc(thread uint4 &st) {
    float c = 2.0f * rand01(st) - 1.0f, ph = 6.2831853f * rand01(st), sn = sqrt(max(0.0f, 1.0f - c * c));
    float cp, sp = sincos(ph, cp);
    return float3(sn * cp, sn * sp, c);
}

kernel void tc_clear(device uint* cnt [[buffer(8)]], constant TParams& T [[buffer(7)]],
                     device const uint* ctr_in [[buffer(16)]], device uint* ctr_out [[buffer(17)]],
                     device uint* args [[buffer(18)]], uint id [[thread_position_in_grid]])
{
    if ((int)id < T.nblocks) cnt[id] = 0;
    if (id == 0) {
        *ctr_out = 0;
        args[0] = (*ctr_in + 255u) / 256u; args[1] = 1; args[2] = 1;
    }
}

// act_in[0..*ctr_in) are the ions that were live last round; those still live are appended to act_out
// (simdgroup-aggregated atomic) and binned. Position j in act_out indexes blk_of/rank_of.
kernel void tc_count(device const ChIon* ions [[buffer(0)]], constant ChParams& P [[buffer(4)]],
                     device atomic_uint* cnt [[buffer(8)]], device uint* blk_of [[buffer(9)]],
                     device uint* rank_of [[buffer(10)]], constant TParams& T [[buffer(7)]],
                     device const uint* act_in [[buffer(14)]], device uint* act_out [[buffer(15)]],
                     device const uint* ctr_in [[buffer(16)]], device atomic_uint* ctr_out [[buffer(17)]],
                     uint id [[thread_position_in_grid]], uint lane [[thread_index_in_simdgroup]])
{
    uint ion_id = 0; bool active = false; uint b = 0;
    if (id < *ctr_in) {
        ion_id = act_in[id];
        int4 pp = ions[ion_id].p;
        if (pp.w == 0) {
            active = true;
            int3 c = int3(tc_gc(pp.xyz, P)) / BR;
            b = (uint)(c.x + T.nbx * (c.y + T.nby * c.z));
        }
    }
    uint pref = simd_prefix_exclusive_sum(active ? 1u : 0u), tot = simd_sum(active ? 1u : 0u);
    uint base = 0;
    if (lane == 0 && tot > 0) base = atomic_fetch_add_explicit(ctr_out, tot, memory_order_relaxed);
    base = simd_broadcast_first(base);
    if (active) {
        uint pos = base + pref;
        act_out[pos] = ion_id; blk_of[pos] = b;
        rank_of[pos] = atomic_fetch_add_explicit(&cnt[b], 1u, memory_order_relaxed);
    }
}

// Single threadgroup of 1024: exclusive scan of ceil(cnt/TILE_TG) over blocks.
kernel void tc_scan(device const uint* cnt [[buffer(8)]], device uint* offs [[buffer(11)]],
                    device uint* tgb [[buffer(13)]], device uint* args [[buffer(18)]],
                    device uint* stats [[buffer(19)]], constant TParams& T [[buffer(7)]],
                    uint lid [[thread_index_in_threadgroup]], uint lane [[thread_index_in_simdgroup]],
                    uint sid [[simdgroup_index_in_threadgroup]], uint nsimd [[simdgroups_per_threadgroup]],
                    uint nthr [[threads_per_threadgroup]])
{
    threadgroup uint sh_tg[32], sh_act[32];
    uint nb = (uint)T.nblocks;
    uint per = (nb + nthr - 1) / nthr;
    uint b0 = min(lid * per, nb), b1 = min(b0 + per, nb);
    uint mysum = 0, myact = 0;
    for (uint b = b0; b < b1; ++b) { uint c = cnt[b]; mysum += (c + TILE_TG - 1) / TILE_TG; myact += c; }
    uint within = simd_prefix_exclusive_sum(mysum);
    uint ssum = simd_sum(mysum), sact = simd_sum(myact);
    if (lane == 0) { sh_tg[sid] = ssum; sh_act[sid] = sact; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    uint so = 0, total_tg = 0, total_act = 0;
    for (uint j = 0; j < nsimd; ++j) { if (j < sid) so += sh_tg[j]; total_tg += sh_tg[j]; total_act += sh_act[j]; }
    uint base = so + within;
    for (uint b = b0; b < b1; ++b) {
        uint c = cnt[b], n = (c + TILE_TG - 1) / TILE_TG;
        offs[b] = base * TILE_TG;
        for (uint j = 0; j < n; ++j) tgb[base + j] = b;
        base += n;
    }
    if (lid == 0) {
        args[4] = (total_act + 255u) / 256u; args[5] = 1; args[6] = 1;     // scatter
        args[8] = total_tg; args[9] = 1; args[10] = 1;                      // push
        stats[0] = total_tg; stats[1] = total_act; stats[2] += total_act; stats[3] += total_tg;
    }
}

kernel void tc_scatter(device const uint* blk_of [[buffer(9)]], device const uint* rank_of [[buffer(10)]],
                       device const uint* offs [[buffer(11)]], device uint* sorted [[buffer(12)]],
                       device const uint* act_out [[buffer(15)]], device const uint* ctr_out [[buffer(17)]],
                       uint id [[thread_position_in_grid]])
{
    if (id >= *ctr_out) return;
    sorted[offs[blk_of[id]] + rank_of[id]] = act_out[id];
}

[[max_total_threads_per_threadgroup(TILE_TG)]]
kernel void tc_push(device ChIon* ions [[buffer(0)]], device const float4* F [[buffer(1)]],
                    constant ChParams& P [[buffer(4)]], device const ChJob* jobs [[buffer(5)]],
                    device const uint* jidb [[buffer(6)]], constant TParams& T [[buffer(7)]],
                    device const uint* cnt [[buffer(8)]], device const uint* sorted [[buffer(12)]],
                    device const uint* offs [[buffer(11)]], device const uint* tgb [[buffer(13)]],
                    device atomic_uint* stats [[buffer(19)]],
                    uint tg [[threadgroup_position_in_grid]], uint lid [[thread_index_in_threadgroup]])
{
    threadgroup float4 sm_rf[BNODES];
    threadgroup float4 sm_dc[BNODES];
    threadgroup float4 sm_ac[BNODES];

    const uint blk = tgb[tg];
    const int bxi = (int)blk % T.nbx, byi = ((int)blk / T.nbx) % T.nby, bzi = (int)blk / (T.nbx * T.nby);
    const int3 cell0 = int3(bxi, byi, bzi) * BR;
    {
        device const float4* src = F + (uint)blk * (BNODES * 3);
        for (uint e = lid; e < BNODES; e += TILE_TG) { sm_rf[e] = src[3 * e]; sm_dc[e] = src[3 * e + 1]; sm_ac[e] = src[3 * e + 2]; }
    }
    const uint slot = tg * TILE_TG + lid;
    const bool valid = (slot - offs[blk]) < cnt[blk];
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (!valid) return;

    const uint id = sorted[slot];
    ChIon ion = ions[id];
    const ChJob J = jobs[jidb[id]];
    int3 p = ion.p.xyz; float3 v = ion.v.xyz; float t = ion.v.w; float4 aux = ion.aux; uint4 rs = ion.rng;
    const float inv_scale = 1.0f / P.pos_scale;
    const float3 lim = float3((P.nx - 1) * P.h, (P.ny - 1) * P.h, (P.nz - 1) * P.h);
    const float3 gmax = float3(P.nx - 1.001f, P.ny - 1.001f, P.nz - 1.001f);
    int status = 0;
    int s = 0;
    for (; s < T.chunk; ++s) {
        float3 x = float3(p) * inv_scale;
        float3 gc = clamp(x * P.inv_h, float3(0.0f), gmax);
        int3 i = int3(gc);
        int3 li = i - cell0;
        if (s == 0) { li = clamp(li, int3(0), int3(BR - 1)); i = li + cell0; }
        else if (any(li < 0) || any(li >= int3(BR))) break;      // park: state untouched
        int b = li.x + BN * (li.y + BN * li.z);
        float a_rf = P.v_rf * sinpi(2.0f * aux.w);
        float a_ac = J.v_ac * sinpi(2.0f * aux.x);
#if CH_NEAREST
        float3 E = a_rf * sm_rf[b].xyz + sm_dc[b].xyz + a_ac * sm_ac[b].xyz;
#else
        float3 tt = gc - float3(i);
        float3 E = a_rf * samp_tri<1>((threadgroup const float4*)sm_rf, b, 0, tt) + samp_tri<1>((threadgroup const float4*)sm_dc, b, 0, tt)
                 + a_ac * samp_tri<1>((threadgroup const float4*)sm_ac, b, 0, tt);
#endif
        // ---- physics: verbatim copy of `channel` ----
        v += (J.qm_kick * P.dt) * E;
        p += int3(rint(v * (P.dt * P.pos_scale)));
        t += P.dt;
        aux.w += P.w_rf; aux.w -= floor(aux.w);
        aux.x += J.w_ac; aux.x -= floor(aux.x);
        if (J.legacy_floor < 0.0f) {
            float bound = length(v) + J.v_bound;
            if (rand01(rs) < J.nsig_dt * bound) {
                float3 vg = J.s_gas * gauss3(rs);
                float3 gv = v - vg; float gl = length(gv);
                if (rand01(rs) * bound < gl) { v = v - J.mr * gv + J.mr * gl * iso_dir(rs); aux.y += 1.0f; }
            }
        } else {
            aux.z -= 1.0f;
            if (aux.z <= 0.0f) {
                float vm = length(v);
                if (vm > J.legacy_floor) {
                    float3 d = iso_dir(rs);
                    float3 vcm = v * (1.0f - J.mr);
                    v = vcm + J.mr * vm * d;
                    aux.y += 1.0f;
                }
                float mfp = -log(max(rand01(rs), 1e-7f)) * J.legacy_mfp_mean;
                aux.z = max(1.0f, floor(mfp / max(length(v), 1e-3f) / P.dt));
            }
        }
        float3 xn = float3(p) * inv_scale;
        if (xn.z < P.z_return) { status = 2; break; }
        if (xn.z > P.z_exit)   { status = 3; break; }
        if (any(xn < 0.0f) || any(xn > lim)) { status = 4; break; }
        if (xn.z > P.crash_z0 && xn.z < P.crash_z1) {
            int3 nn = int3(rint(xn * P.inv_h));
            int3 ln = nn - cell0;
            float mk;
            if (all(ln >= 0) && all(ln < int3(BN))) mk = sm_rf[ln.x + BN * (ln.y + BN * ln.z)].w;
            else mk = F[brick_node(nn, T.nbx, T.nby)].w;
            if (mk > 0.5f) { status = 1; break; }
        }
        if (t > P.t_max) { status = 5; break; }
    }
    ion.p = int4(p, status); ion.v = float4(v, t); ion.aux = aux; ion.rng = rs;
    ions[id] = ion;
#ifdef CH_STATS
    uint ss = simd_sum((uint)s);
    if (simd_is_first()) atomic_fetch_add_explicit(&stats[4], ss, memory_order_relaxed);
#endif
}

// ---------------------------------------------------------------------------------------------------------------
// Tail: compaction + persistent workers.
// ---------------------------------------------------------------------------------------------------------------

// Compact the live ions of act_in[0..*ctr_in) into act_out (dense lanes for the workers).
kernel void tc_filter(device const ChIon* ions [[buffer(0)]], device const uint* act_in [[buffer(14)]],
                      device uint* act_out [[buffer(15)]], device const uint* ctr_in [[buffer(16)]],
                      device atomic_uint* ctr_out [[buffer(17)]], uint id [[thread_position_in_grid]],
                      uint lane [[thread_index_in_simdgroup]])
{
    uint ion_id = 0; bool active = false;
    if (id < *ctr_in) { ion_id = act_in[id]; active = ions[ion_id].p.w == 0; }
    uint pref = simd_prefix_exclusive_sum(active ? 1u : 0u), tot = simd_sum(active ? 1u : 0u);
    uint base = 0;
    if (lane == 0 && tot > 0) base = atomic_fetch_add_explicit(ctr_out, tot, memory_order_relaxed);
    base = simd_broadcast_first(base);
    if (active) act_out[base + pref] = ion_id;
}
kernel void tc_zero_args(device uint* ctr_out [[buffer(17)]], device const uint* ctr_in [[buffer(16)]], device uint* args [[buffer(18)]],
                         device uint* head [[buffer(20)]], uint id [[thread_position_in_grid]])
{
    if (id == 0) { *ctr_out = 0; *head = 0; args[0] = (*ctr_in + 255u) / 256u; args[1] = 1; args[2] = 1; }
}
// Worker launch geometry from the live count: at most TAIL_WMAX workers; when few ions are alive, spread them over
// simdgroups with only `L` active lanes each (less divergence and less contention => lower per-step latency).
kernel void tc_size_args(device const uint* ctr_out [[buffer(17)]], device uint* args [[buffer(18)]], device uint* cfg [[buffer(21)]],
                         uint id [[thread_position_in_grid]])
{
    if (id == 0) {
        uint count = *ctr_out;
        uint L = 1;
        while (L < 32u && count > TAIL_LSPAN * L) L *= 2u;
        uint workers = min(count, (uint)TAIL_WMAX);
        uint groups = (workers + L - 1u) / L;
        uint threads = groups * 32u;
        args[4] = (threads + TAIL_TG - 1u) / TAIL_TG; args[5] = 1; args[6] = 1;
        cfg[0] = L;
    }
}

// Persistent-worker per-ion kernel (direct, untiled, brick-layout lookups). Each worker repeatedly takes the next live ion from
// the compacted list and advances it by up to P.nsteps steps, so lanes never idle while work remains; ions that outlive the
// chunk are picked up by the next launch.
//
// Latency hiding (arithmetic unchanged): the field data, the conductor-mask node and the RF/AC drive of the NEXT step depend only
// on the new position and phases, never on the collision outcome, so their loads / sinpi are issued right after the position
// update, before the collision and termination code (which only touch v, rng, aux.y/z and the status).
#ifdef TAIL_MAXT
[[max_total_threads_per_threadgroup(TAIL_MAXT)]]
#endif
kernel void tc_tail(device ChIon* ions [[buffer(0)]], device const float4* F [[buffer(1)]],
                    constant ChParams& P [[buffer(4)]], device const ChJob* jobs [[buffer(5)]],
                    device const uint* jidb [[buffer(6)]], device const uint* act [[buffer(14)]],
                    device const uint* ctr [[buffer(17)]], device atomic_uint* head [[buffer(20)]],
                    device const uint* cfg [[buffer(21)]], uint lane [[thread_index_in_simdgroup]])
{
    if (lane >= cfg[0]) return;
    const uint count = *ctr;
    const int nbx = NBX, nby = NBY;
    const float inv_scale = 1.0f / P.pos_scale;
    const float3 lim = float3((P.nx - 1) * P.h, (P.ny - 1) * P.h, (P.nz - 1) * P.h);
    const float3 gmax = float3(P.nx - 1.001f, P.ny - 1.001f, P.nz - 1.001f);
    const int3 nmax = int3(P.nx - 1, P.ny - 1, P.nz - 1);
    for (;;) {
        const uint idx = atomic_fetch_add_explicit(head, 1u, memory_order_relaxed);
        if (idx >= count) return;
        const uint id = act[idx];
        ChIon ion = ions[id];
        if (ion.p.w != 0) continue;
        const float jq = jobs[jidb[id]].qm_kick, jv_ac = jobs[jidb[id]].v_ac, jw_ac = jobs[jidb[id]].w_ac, jnsig = jobs[jidb[id]].nsig_dt, jsgas = jobs[jidb[id]].s_gas,
                    jvb = jobs[jidb[id]].v_bound, jmr = jobs[jidb[id]].mr, jfloor = jobs[jidb[id]].legacy_floor, jmfp = jobs[jidb[id]].legacy_mfp_mean;
        int3 p = ion.p.xyz; float3 v = ion.v.xyz; float t = ion.v.w; float4 aux = ion.aux; uint4 rs = ion.rng;
        int status = 0;
        // field data for the current position
        float3 gc, tt; int3 i; int b;
        {
            float3 x = float3(p) * inv_scale;
            gc = clamp(x * P.inv_h, float3(0.0f), gmax); i = int3(gc); b = brick_node(i, nbx, nby);
        }
#if CH_NEAREST
        float3 f_rf = F[b].xyz, f_dc = F[b + 1].xyz, f_ac = F[b + 2].xyz;
#else
        CellCache cc; cell_load(cc, F, b, 3, 1, i); tt = gc - float3(i);
#endif
        float a_rf = P.v_rf * sinpi(2.0f * aux.w);
        float a_ac = jv_ac * sinpi(2.0f * aux.x);
        for (int s = 0; s < P.nsteps; ++s) {
#if CH_NEAREST
            float3 E = a_rf * f_rf + f_dc + a_ac * f_ac;
#else
#ifdef AC_FULL
            float3 E = a_rf * cell_mix(cc.c[0], tt) + cell_mix(cc.c[1], tt) + a_ac * cell_mix(cc.c[2], tt);
#else
            float3 E = a_rf * mix8(cc.r0, cc.r1, cc.r2, cc.r3, cc.r4, cc.r5, cc.r6, cc.r7, tt) + mix8(cc.d0, cc.d1, cc.d2, cc.d3, cc.d4, cc.d5, cc.d6, cc.d7, tt) + a_ac * cc.ac;
#endif
#endif
            v += (jq * P.dt) * E;
            p += int3(rint(v * (P.dt * P.pos_scale)));
            t += P.dt;
            aux.w += P.w_rf; aux.w -= floor(aux.w);
            aux.x += jw_ac; aux.x -= floor(aux.x);
            // ---- everything the next step / the termination test needs that depends only on p and the phases ----
            const float3 xn = float3(p) * inv_scale;
            const float3 g = xn * P.inv_h;
            gc = clamp(g, float3(0.0f), gmax); i = int3(gc); b = brick_node(i, nbx, nby);
#if CH_NEAREST
            f_rf = F[b].xyz; f_dc = F[b + 1].xyz; f_ac = F[b + 2].xyz;
#else
            tt = gc - float3(i);
            if (any(i != cc.cell)) cell_load(cc, F, b, 3, 1, i);
#endif
            float mk = 0.0f;                                      // conductor mask at the nearest node (clamped address; used only if in range)
            if (xn.z > P.crash_z0 && xn.z < P.crash_z1) mk = F[brick_node(clamp(int3(rint(g)), int3(0), nmax), nbx, nby)].w;
            a_rf = P.v_rf * sinpi(2.0f * aux.w);
            a_ac = jv_ac * sinpi(2.0f * aux.x);
            // ---- collisions (as in `channel`) ----
            if (jfloor < 0.0f) {
                float bound = length(v) + jvb;
                if (rand01(rs) < jnsig * bound) {
                    float3 vg = jsgas * gauss3_sc(rs);
                    float3 gv = v - vg; float gl = length(gv);
                    if (rand01(rs) * bound < gl) { v = v - jmr * gv + jmr * gl * iso_dir_sc(rs); aux.y += 1.0f; }
                }
            } else {
                aux.z -= 1.0f;
                if (aux.z <= 0.0f) {
                    float vm = length(v);
                    if (vm > jfloor) {
                        float3 d = iso_dir_sc(rs);
                        float3 vcm = v * (1.0f - jmr);
                        v = vcm + jmr * vm * d;
                        aux.y += 1.0f;
                    }
                    float mfp = -log(max(rand01(rs), 1e-7f)) * jmfp;
                    aux.z = max(1.0f, floor(mfp / max(length(v), 1e-3f) / P.dt));
                }
            }
            if (xn.z < P.z_return) { status = 2; break; }
            if (xn.z > P.z_exit)   { status = 3; break; }
            if (any(xn < 0.0f) || any(xn > lim)) { status = 4; break; }
            if (xn.z > P.crash_z0 && xn.z < P.crash_z1 && mk > 0.5f) { status = 1; break; }
            if (t > P.t_max) { status = 5; break; }
        }
        ion.p = int4(p, status); ion.v = float4(v, t); ion.aux = aux; ion.rng = rs;
        ions[id] = ion;
    }
}

// ---------------------------------------------------------------------------------------------------------------
// Legacy-collision worker with collision batching (the legacy model collides every ~4 steps, so in lock-step a simdgroup would run
// the expensive collision code (sincos, log, sqrt, divides) on nearly every step with ~25% of its lanes active).
// Here a lane that reaches its collision step waits (`pending`) while the others keep stepping, and the collision code runs once
// enough lanes are waiting. Each ion still executes exactly the sequence of the reference kernel (step, collision, termination
// test; no ion step happens between its collision-step and its collision), so states are bitwise identical.
// ---------------------------------------------------------------------------------------------------------------
kernel void tc_tail_leg(device ChIon* ions [[buffer(0)]], device const float4* F [[buffer(1)]],
                        constant ChParams& P [[buffer(4)]], device const ChJob* jobs [[buffer(5)]],
                        device const uint* jidb [[buffer(6)]], device const uint* act [[buffer(14)]],
                        device const uint* ctr [[buffer(17)]], device atomic_uint* head [[buffer(20)]],
                        device const uint* cfg [[buffer(21)]], uint lane [[thread_index_in_simdgroup]])
{
    if (lane >= cfg[0]) return;
    const uint count = *ctr;
    const int nbx = NBX, nby = NBY;
    const float inv_scale = 1.0f / P.pos_scale;
    const float3 lim = float3((P.nx - 1) * P.h, (P.ny - 1) * P.h, (P.nz - 1) * P.h);
    const float3 gmax = float3(P.nx - 1.001f, P.ny - 1.001f, P.nz - 1.001f);
    const int3 nmax = int3(P.nx - 1, P.ny - 1, P.nz - 1);
    bool have = false, pending = false, finishing = false, drained = false;
    uint id = 0; ChIon ion; ChJob J = jobs[0];
    int3 p = int3(0); float3 v = float3(0.0f); float t = 0.0f; float4 aux = float4(0.0f); uint4 rs = uint4(1u);
    int status = 0, s = 0;
    float3 gc = float3(0.0f), tt = float3(0.0f); int3 i = int3(0); int b = 0;
    float3 f_rf = float3(0.0f), f_dc = float3(0.0f), f_ac = float3(0.0f);
    float a_rf = 0.0f, a_ac = 0.0f;
    for (;;) {
        // (1) acquire an ion
        if (!have && !drained) {
            const uint idx = atomic_fetch_add_explicit(head, 1u, memory_order_relaxed);
            if (idx >= count) { drained = true; }
            else {
                id = act[idx]; ion = ions[id];
                if (ion.p.w == 0) {
                    have = true; pending = false; finishing = false; status = 0; s = 0;
                    J = jobs[jidb[id]];
                    p = ion.p.xyz; v = ion.v.xyz; t = ion.v.w; aux = ion.aux; rs = ion.rng;
                    float3 x = float3(p) * inv_scale;
                    gc = clamp(x * P.inv_h, float3(0.0f), gmax); i = int3(gc); b = brick_node(i, nbx, nby);
                    f_rf = F[b].xyz; f_dc = F[b + 1].xyz; f_ac = F[b + 2].xyz;
                    a_rf = P.v_rf * sinpi(2.0f * aux.w);
                    a_ac = J.v_ac * sinpi(2.0f * aux.x);
                }
            }
        }
        if (!simd_any(have || !drained)) break;
        // (2) one step for every lane that holds an ion and is not waiting for its collision
        if (have && !pending && !finishing) {
            float3 E = a_rf * f_rf + f_dc + a_ac * f_ac;
            v += (J.qm_kick * P.dt) * E;
            p += int3(rint(v * (P.dt * P.pos_scale)));
            t += P.dt;
            aux.w += P.w_rf; aux.w -= floor(aux.w);
            aux.x += J.w_ac; aux.x -= floor(aux.x);
            const float3 xn = float3(p) * inv_scale;
            const float3 g = xn * P.inv_h;
            gc = clamp(g, float3(0.0f), gmax); i = int3(gc); b = brick_node(i, nbx, nby);
            f_rf = F[b].xyz; f_dc = F[b + 1].xyz; f_ac = F[b + 2].xyz;
            float mk = 0.0f;
            if (xn.z > P.crash_z0 && xn.z < P.crash_z1) mk = F[brick_node(clamp(int3(rint(g)), int3(0), nmax), nbx, nby)].w;
            a_rf = P.v_rf * sinpi(2.0f * aux.w);
            a_ac = J.v_ac * sinpi(2.0f * aux.x);
            aux.z -= 1.0f;
            pending = aux.z <= 0.0f;                              // collision of this step is applied in phase (3)
            // termination test of this step (depends on p and t only, so it can precede the deferred collision)
            if (xn.z < P.z_return) status = 2;
            else if (xn.z > P.z_exit) status = 3;
            else if (any(xn < 0.0f) || any(xn > lim)) status = 4;
            else if (xn.z > P.crash_z0 && xn.z < P.crash_z1 && mk > 0.5f) status = 1;
            else if (t > P.t_max) status = 5;
            ++s;
            finishing = (status != 0) || (s >= P.nsteps);
        }
        // (3) collisions, batched: run the collision code once at least half of the lanes holding ions are waiting
        {
            const uint npend = simd_sum(pending ? 1u : 0u);
            const uint nhave = simd_sum(have ? 1u : 0u);
            const uint nrun = simd_sum((have && !pending && !finishing) ? 1u : 0u);
            if (npend > 0 && (2u * npend >= nhave || nrun == 0)) {
                if (pending) {
                    float vm = length(v);
                    if (vm > J.legacy_floor) {
                        float3 d = iso_dir_sc(rs);
                        float3 vcm = v * (1.0f - J.mr);
                        v = vcm + J.mr * vm * d;
                        aux.y += 1.0f;
                    }
                    float mfp = -log(max(rand01(rs), 1e-7f)) * J.legacy_mfp_mean;
                    aux.z = max(1.0f, floor(mfp / max(length(v), 1e-3f) / P.dt));
                    pending = false;
                }
            }
        }
        // (4) write back finished ions
        if (have && finishing && !pending) {
            ion.p = int4(p, status); ion.v = float4(v, t); ion.aux = aux; ion.rng = rs;
            ions[id] = ion;
            have = false;
        }
    }
}
