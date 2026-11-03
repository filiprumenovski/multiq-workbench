// Ejection-channel kernel: per-ion time/phase, outcome classification, legacy switches.
#include <metal_stdlib>
using namespace metal;

// Batch-wide parameters (shared by every ion in a launch).
struct ChParams {
    float dt;
    float v_rf, w_rf;           // w_* hold f*dt (cycles per step)
    float inv_h, h;
    float pos_scale;
    float z_return, z_exit;     // m
    float t_max;                // s
    float crash_z0, crash_z1;   // conductor-hit test only for crash_z0 < z < crash_z1 (legacy: rod region only)
    int nx, ny, nz;
    int nsteps;
    int nearest;                // 1 = legacy nearest-node field sampling
    int legacy;                 // 1 = legacy stationary-gas collisions in every job (selects the batched-collision worker kernel)
};

// Per-job (per-species / per-AC-drive) parameters; each ion carries a job id (parallel buffer `jid`).
struct ChJob {
    float qm_kick;              // q/m (x0.5 in legacy half-force mode)
    float v_ac, w_ac;           // AC amplitude and cycles per step
    float nsig_dt;              // n*sigma*dt (correct mode)
    float s_gas, v_bound, mr;   // Maxwell gas std, null-collision bound add-on, M/(m+M)
    float legacy_floor;         // legacy: no collision below this speed [m/s]; <0 => correct thermal collisions
    float legacy_mfp_mean;      // legacy: mean free path [m] (collision every mfp/|v|)
    float pad0, pad1, pad2;
};

struct ChIon {
    int4 p;        // fixed-point pos; w = status: 0 alive, 1 crashed, 2 returned, 3 transmitted, 4 left box, 5 timeout
    float4 v;      // velocity; w = own time [s]
    float4 aux;    // x = AC phase [cycles, accumulating in 0..1), y = collisions, z = legacy steps-to-next-collision, w = RF phase [cycles]
    uint4 rng;
};

static inline uint rotl(uint x, int k) { return (x << k) | (x >> (32 - k)); }
static inline uint next_u32(thread uint4 &s) {
    uint r = rotl(s.y * 5u, 7) * 9u; uint t = s.y << 9;
    s.z ^= s.x; s.w ^= s.y; s.y ^= s.z; s.x ^= s.w; s.z ^= t; s.w = rotl(s.w, 11);
    return r;
}
static inline float rand01(thread uint4 &s) { return (float)(next_u32(s) >> 8) * (1.0f / 16777216.0f); }
static inline float3 gauss3(thread uint4 &s) {
    float u1 = max(rand01(s), 1e-7f), u2 = rand01(s), u3 = max(rand01(s), 1e-7f), u4 = rand01(s);
    float r1 = sqrt(-2.0f * log(u1)), r2 = sqrt(-2.0f * log(u3));
    return float3(r1 * cos(6.2831853f * u2), r1 * sin(6.2831853f * u2), r2 * cos(6.2831853f * u4));
}
static inline float3 iso_dir(thread uint4 &s) {
    float c = 2.0f * rand01(s) - 1.0f, ph = 6.2831853f * rand01(s), sn = sqrt(max(0.0f, 1.0f - c * c));
    return float3(sn * cos(ph), sn * sin(ph), c);
}

struct Tri { int b; float3 t; };
static inline Tri locate(float3 g, int nx, int ny, int nz, int nearest) {
    float3 gc = clamp(g, float3(0.0f), float3(nx - 1.001f, ny - 1.001f, nz - 1.001f));
    int3 i = int3(gc);
    Tri r; r.b = i.x + nx * (i.y + ny * i.z);
    r.t = nearest ? float3(0.0f) : gc - float3(i);
    return r;
}
static inline float3 interp(device const float4* F, Tri c, int nx, int ny) {
    int sy = nx, sz = nx * ny, b = c.b; float3 t = c.t;
    float3 c00 = mix(F[b].xyz, F[b + 1].xyz, t.x),           c10 = mix(F[b + sy].xyz, F[b + 1 + sy].xyz, t.x);
    float3 c01 = mix(F[b + sz].xyz, F[b + 1 + sz].xyz, t.x), c11 = mix(F[b + sy + sz].xyz, F[b + 1 + sy + sz].xyz, t.x);
    return mix(mix(c00, c10, t.y), mix(c01, c11, t.y), t.z);
}

kernel void channel(device ChIon* ions [[buffer(0)]], device const float4* Frf [[buffer(1)]],
                    device const float4* Fdc [[buffer(2)]], device const float4* Fac [[buffer(3)]],
                    constant ChParams& P [[buffer(4)]], device const ChJob* jobs [[buffer(5)]],
                    device const uint* jid [[buffer(6)]], uint id [[thread_position_in_grid]])
{
    ChIon ion = ions[id];
    if (ion.p.w != 0) return;
    const ChJob J = jobs[jid[id]];
    int3 p = ion.p.xyz; float3 v = ion.v.xyz; float t = ion.v.w; float4 aux = ion.aux; uint4 rs = ion.rng;
    const float inv_scale = 1.0f / P.pos_scale;
    const float3 lim = float3((P.nx - 1) * P.h, (P.ny - 1) * P.h, (P.nz - 1) * P.h);
    int status = 0;
    for (int s = 0; s < P.nsteps; ++s) {
        float3 x = float3(p) * inv_scale;
        Tri c = locate(x * P.inv_h, P.nx, P.ny, P.nz, P.nearest);
        float a_rf = P.v_rf * sinpi(2.0f * aux.w);       // phases kept as cycle accumulators in [0,1):
        float a_ac = J.v_ac * sinpi(2.0f * aux.x);       // no float32 time-accumulation phase drift
        float3 E = a_rf * interp(Frf, c, P.nx, P.ny) + interp(Fdc, c, P.nx, P.ny) + a_ac * interp(Fac, c, P.nx, P.ny);
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
            if (aux.z <= 0.0f) {                       // legacy: collide with stationary gas if fast enough
                float vm = length(v);
                if (vm > J.legacy_floor) {
                    float3 d = iso_dir(rs);             // isotropic CM scattering off gas at rest
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
        int3 nn = int3(rint(xn * P.inv_h));
        if (xn.z > P.crash_z0 && xn.z < P.crash_z1 && Frf[nn.x + P.nx * (nn.y + P.ny * nn.z)].w > 0.5f) { status = 1; break; }
        if (t > P.t_max) { status = 5; break; }
    }
    ion.p = int4(p, status); ion.v = float4(v, t); ion.aux = aux; ion.rng = rs;
    ions[id] = ion;
}
