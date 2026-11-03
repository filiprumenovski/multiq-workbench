//! Small parallel helpers shared by the kernels.
//!
//! The multigrid kernels split work into contiguous "lines" (runs of `nx`
//! nodes) and hand contiguous ranges of lines to rayon tasks.  Writes from
//! different tasks always go to disjoint memory (or, for the red-black
//! smoother, to nodes of one colour while only nodes of the other colour are
//! read), so the kernels share raw pointers through [`SyncPtr`].

use rayon::prelude::*;

/// Raw mutable pointer that may be shared between rayon tasks.
///
/// # Safety contract
/// The creator guarantees that concurrent users never write to the same
/// element and never read an element another task writes.
#[derive(Clone, Copy)]
pub(crate) struct SyncPtr<T>(*mut T);
unsafe impl<T> Send for SyncPtr<T> {}
unsafe impl<T> Sync for SyncPtr<T> {}

impl<T> SyncPtr<T> {
    #[inline]
    pub fn new(s: &mut [T]) -> Self {
        SyncPtr(s.as_mut_ptr())
    }
    /// Accessor (a method, so closures capture the whole wrapper, not the raw field).
    #[inline(always)]
    pub fn get(self) -> *mut T {
        self.0
    }
}

/// Minimum amount of work (in nodes) below which we stay on one thread.
const SERIAL_WORK: usize = 1 << 15;

/// Number of chunks for `n` items of `cost` work each.
fn nchunks(n: usize, cost: usize) -> usize {
    let threads = rayon::current_num_threads();
    if threads <= 1 || n.saturating_mul(cost) < SERIAL_WORK {
        return 1;
    }
    // Several chunks per thread so rayon can balance irregular (masked) work.
    (threads * 6).min(n).max(1)
}

/// Run `f(start, end)` over a partition of `0..n`; `cost` is the approximate
/// work per item and only steers the serial/parallel decision.
pub(crate) fn par_ranges<F>(n: usize, cost: usize, f: F)
where
    F: Fn(usize, usize) + Sync + Send,
{
    if n == 0 {
        return;
    }
    let nc = nchunks(n, cost);
    if nc == 1 {
        f(0, n);
        return;
    }
    let chunk = n.div_ceil(nc);
    (0..nc).into_par_iter().for_each(|c| {
        let a = c * chunk;
        let b = ((c + 1) * chunk).min(n);
        if a < b {
            f(a, b);
        }
    });
}

/// Like [`par_ranges`] but each range returns an `f64`; results are max-reduced.
pub(crate) fn par_max<F>(n: usize, cost: usize, f: F) -> f64
where
    F: Fn(usize, usize) -> f64 + Sync + Send,
{
    if n == 0 {
        return 0.0;
    }
    let nc = nchunks(n, cost);
    if nc == 1 {
        return f(0, n);
    }
    let chunk = n.div_ceil(nc);
    (0..nc)
        .into_par_iter()
        .map(|c| {
            let a = c * chunk;
            let b = ((c + 1) * chunk).min(n);
            if a < b { f(a, b) } else { 0.0 }
        })
        .reduce(|| 0.0, f64::max)
}
