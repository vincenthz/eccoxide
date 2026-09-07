//! Batched Bernstein-Yang ("safegcd") constant-time modular inversion.
//!
//! The inversion is a sequence of *divsteps* on a state `(delta, f, g)` paired
//! with two accumulators `(v, r)`:
//!
//! ```text
//! if delta > 0 and g is odd:  (delta, f, g) <- (1 - delta, g, (g - f) / 2)
//! otherwise:                  (delta, f, g) <- (1 + delta, f, (g + (g & 1) * f) / 2)
//! ```
//!
//! started from `(1, m, a)`, `(v, r) = (0, 1)`, and iterated a number of times
//! fixed by the bit length of the modulus. Every step preserves
//!
//! ```text
//! v * a = f * 2^k  (mod m)      and      r * a = g * 2^k  (mod m)
//! ```
//!
//! (`k` being the number of steps performed), and after enough steps `g` is
//! zero and `f` is `±1`, so `1/a = sign(f) * v * 2^-k`.
//!
//! fiat-crypto emits this as a single-step divstep primitive. The point
//! of this module is that a *batch* of [`BATCH`] consecutive divsteps only
//! inspects the low bit of `f` and `g` at each step, so the whole batch can be
//! run on the low limbs alone, in registers, accumulating a 2x2 integer
//! transition matrix. That matrix is then applied to the full-width values
//! once per batch:
//!
//! ```text
//! 2^BATCH * f' = ff * f + fg * g          v' = (ff * v + fg * r) / 2^BATCH
//! 2^BATCH * g' = gf * f + gg * g          r' = (gf * v + gg * r) / 2^BATCH
//! ```
//!
//! For a 256-bit modulus, it means turning ~741 full-width passes into 12.
//!
//! The accumulators divide by `2^BATCH` per batch (a Montgomery-style
//! reduction, [`update_acc`]) rather than multiplying up, which keeps them
//! reduced and makes the total division `2^(BATCH * batches)` exactly cancel
//! the `2^k` in the invariant. No final correction constant is therefore
//! needed beyond the sign of `f`.
//!
//! Everything here is constant time with respect to the value being inverted:
//! the batch count depends only on the modulus, and every data-dependent
//! choice is a mask, never a branch.

use crate::mp::ct::{Choice, CtZero};

/// Number of divsteps performed per batch.
///
/// A batch runs on the low 64-bit limbs of `f` and `g`; each step shifts them
/// right by one and pulls in a zero rather than the true next bit, so after
/// `k` steps only the low `64 - k` bits are correct. Step `k` reads bit 0, so
/// any `BATCH <= 63` is sound. 62 is the usual choice: it leaves two bits of
/// slack and bounds the matrix entries by `2^62`, which keeps them (and their
/// doubling) inside `i64`.
const BATCH: u32 = 62;

const MASK: u64 = (1u64 << BATCH) - 1;

/// The 2x2 integer transition matrix of a batch of divsteps.
///
/// Rows are `(ff, fg)` for `f` and `(gf, gg)` for `g`. Bernstein-Yang bounds
/// `|ff| + |fg|` and `|gf| + |gg|` by `2^BATCH`.
#[derive(Clone, Copy)]
struct TransitionMatrix {
    ff: i64,
    fg: i64,
    gf: i64,
    gg: i64,
}

/// Run [`BATCH`] divsteps on the low limbs, returning the new `eta` and the
/// transition matrix. Branch-free in `f` and `g`.
///
/// `eta` is `-delta`, which makes "is a swap possible" (`delta > 0`) just the
/// sign bit. The swap itself never materialises as a select: negating `f` and
/// the `f`-row and folding them into `g` and the `g`-row, then folding the
/// *updated* `g` back into `f`, produces both branches of the divstep. With
/// `c1` the swap mask and `c2` the "g is odd" mask:
///
/// ```text
/// swap:      g += -f  =>  g - f ;  f += (g - f)  =>  g     ; then g >>= 1
/// g odd:     g +=  f  =>  g + f ;  f unchanged             ; then g >>= 1
/// g even:    g, f unchanged                                ; then g >>= 1
/// ```
fn divsteps(mut eta: i64, mut f: u64, mut g: u64) -> (i64, TransitionMatrix) {
    let (mut u, mut v, mut q, mut r) = (1i64, 0i64, 0i64, 1i64);

    let mut i = 0;
    while i < BATCH {
        let c1 = eta >> 63; // all ones when a swap is possible
        let c2 = (g & 1).wrapping_neg() as i64; // all ones when g is odd

        // ±f and +/-(f-row), folded into g and the g-row
        let x = ((f as i64) ^ c1).wrapping_sub(c1);
        let y = (u ^ c1).wrapping_sub(c1);
        let z = (v ^ c1).wrapping_sub(c1);
        g = (g as i64).wrapping_add(x & c2) as u64;
        q = q.wrapping_add(y & c2);
        r = r.wrapping_add(z & c2);

        // a swap needs g odd too; fold the updated g-row back into the f-row
        let c1 = c1 & c2;
        eta = (eta ^ c1).wrapping_sub(c1.wrapping_add(1));
        f = f.wrapping_add(g & (c1 as u64));
        u = u.wrapping_add(q & c1);
        v = v.wrapping_add(r & c1);

        g >>= 1;
        u <<= 1;
        v <<= 1;

        i += 1;
    }

    (
        eta,
        TransitionMatrix {
            ff: u,
            fg: v,
            gf: q,
            gg: r,
        },
    )
}

/// `-m^-1 mod 2^BATCH`, from the low limb of the (odd) modulus.
///
/// Newton iteration doubles the number of correct bits each round, so six
/// rounds from a seed correct modulo 2 reach modulo `2^64`.
const fn minv(m0: u64) -> u64 {
    let mut inv: u64 = 1;
    let mut i = 0;
    while i < 6 {
        inv = inv.wrapping_mul(2u64.wrapping_sub(m0.wrapping_mul(inv)));
        i += 1;
    }
    inv.wrapping_neg() & MASK
}

/// `(a * x + b * y) >> BATCH`, truncated to `S` limbs and read as two
/// complement throughout.
///
/// Multiplication modulo `2^(64*S)` does not depend on how the operands are
/// signed, so the limbs of `x` and `y` go in unsigned and only the multipliers
/// carry a sign; the truncation to `S` limbs then yields the right two's
/// complement result, given that the true value fits (Bernstein-Yang keeps
/// `|f|` and `|g|` bounded by the modulus, so one spare limb is enough).
///
/// `|a| + |b| <= 2^BATCH` bounds every accumulator step by
/// `2^BATCH * (2^64 - 1) + 2^63 < 2^127`, so `i128` does not overflow.
fn lincomb_shift<const S: usize>(a: i64, x: &[u64; S], b: i64, y: &[u64; S]) -> [u64; S] {
    let mut prod = [0u64; S];
    let mut acc: i128 = 0;
    let mut i = 0;
    while i < S {
        acc += (a as i128) * (x[i] as i128) + (b as i128) * (y[i] as i128);
        prod[i] = acc as u64;
        acc >>= 64;
        i += 1;
    }

    let mut out = [0u64; S];
    let mut i = 0;
    while i + 1 < S {
        out[i] = (prod[i] >> BATCH) | (prod[i + 1] << (64 - BATCH));
        i += 1;
    }
    out[S - 1] = ((prod[S - 1] as i64) >> BATCH) as u64;
    out
}

/// `acc += s * x`, with `acc` one limb wider than `x`.
fn mul_add_into<const L: usize, const S: usize>(acc: &mut [u64; S], s: u64, x: &[u64; L]) {
    let mut carry: u128 = 0;
    let mut i = 0;
    while i < L {
        let t = (x[i] as u128) * (s as u128) + (acc[i] as u128) + carry;
        acc[i] = t as u64;
        carry = t >> 64;
        i += 1;
    }
    let mut i = L;
    while i < S {
        let t = (acc[i] as u128) + carry;
        acc[i] = t as u64;
        carry = t >> 64;
        i += 1;
    }
}

/// `(|s|, if s < 0 { m - x } else { x })`.
///
/// Folding the sign of `s` into the operand keeps the products unsigned. `x`
/// is allowed to be `m` on the way out (when `x` was zero), which the callers'
/// bounds account for.
fn split_sign<const L: usize>(s: i64, x: &[u64; L], m: &[u64; L]) -> (u64, [u64; L]) {
    let neg = ((s >> 63) as u64) & 1;
    let mask = 0u64.wrapping_sub(neg);
    let abs = ((s ^ mask as i64).wrapping_sub(mask as i64)) as u64;

    let mut out = [0u64; L];
    let mut borrow: u128 = 0;
    let mut i = 0;
    while i < L {
        let t = (m[i] as u128)
            .wrapping_sub(x[i] as u128)
            .wrapping_sub(borrow);
        out[i] = (t as u64 & mask) | (x[i] & !mask);
        borrow = (t >> 64) & 1;
        i += 1;
    }
    (abs, out)
}

/// `(p * v + q * r) / 2^BATCH  (mod m)`, for `v` and `r` in `[0, m]`.
///
/// The division is a Montgomery style reduction: a multiple of `m` is added to
/// clear the low `BATCH` bits before shifting. With `|p| + |q| <= 2^BATCH` the
/// un-reduced sum is below `2^BATCH * m`, the sum plus the multiple of `m` is
/// below `2^(BATCH+1) * m` (so it still fits in `S = L + 1` limbs), and the
/// shifted result is below `2 * m`, which one conditional subtraction reduces.
fn update_acc<const L: usize, const S: usize>(
    m: &[u64; L],
    m_inv: u64,
    p: i64,
    v: &[u64; L],
    q: i64,
    r: &[u64; L],
) -> [u64; L] {
    let (pa, va) = split_sign(p, v, m);
    let (qa, ra) = split_sign(q, r, m);

    let mut acc = [0u64; S];
    mul_add_into(&mut acc, pa, &va);
    mul_add_into(&mut acc, qa, &ra);

    let k = acc[0].wrapping_mul(m_inv) & MASK;
    mul_add_into(&mut acc, k, m);
    debug_assert_eq!(acc[0] & MASK, 0);

    // acc >>= BATCH (logical: acc is non-negative), then reduce once
    let mut sh = [0u64; S];
    let mut i = 0;
    while i + 1 < S {
        sh[i] = (acc[i] >> BATCH) | (acc[i + 1] << (64 - BATCH));
        i += 1;
    }
    sh[S - 1] = acc[S - 1] >> BATCH;

    // t = sh - m, over the full S limbs so the borrow decides the reduction
    let mut t = [0u64; S];
    let mut borrow: u128 = 0;
    let mut i = 0;
    while i < S {
        let mi = if i < L { m[i] } else { 0 };
        let d = (sh[i] as u128)
            .wrapping_sub(mi as u128)
            .wrapping_sub(borrow);
        t[i] = d as u64;
        borrow = (d >> 64) & 1;
        i += 1;
    }
    let mask = 0u64.wrapping_sub(borrow as u64); // all ones when sh < m
    let mut out = [0u64; L];
    let mut i = 0;
    while i < L {
        out[i] = (t[i] & !mask) | (sh[i] & mask);
        i += 1;
    }
    out
}

/// The Bernstein-Yang step count for a modulus given in saturated form.
///
/// It depends only on the bit length of the public modulus, so this
/// leaks nothing about the value being inverted.
pub(crate) fn iterations<const S: usize>(msat: &[u64; S]) -> usize {
    let mut bits = 0usize;
    let mut i = S;
    while i > 0 {
        i -= 1;
        if msat[i] != 0 {
            bits = i * 64 + (64 - msat[i].leading_zeros() as usize);
            break;
        }
    }
    (49 * bits + if bits < 46 { 80 } else { 57 }) / 17
}

/// Modular inverse of `a` modulo the modulus given in saturated form.
///
/// * `msat` is the modulus, little-endian, in `S = L + 1` limbs (as produced by
///   the fiat-crypto `msat` routine): one limb wider than a field element so
///   that `f` and `g` have room for their sign.
/// * `a` is the value to invert, little-endian, reduced, and *not* in the
///   Montgomery domain.
/// * `one` seeds the accumulator with the representation of `1`: pass the
///   Montgomery one to get the result as a Montgomery representative, since
///   the accumulators only ever undergo mod-`m` linear maps.
/// * `iterations` is the Bernstein-Yang step count for this modulus. It is
///   rounded up to a whole number of batches; the extra divsteps are harmless,
///   as once `g` reaches zero the state is stationary apart from the
///   accumulator, whose per-batch division tracks them exactly.
///
/// Returns the accumulator and whether `f` ended up negative, i.e. whether the
/// caller must negate to obtain the inverse.
pub(crate) fn inverse<const L: usize, const S: usize>(
    msat: &[u64; S],
    a: &[u64; L],
    one: &[u64; L],
    iterations: usize,
) -> ([u64; L], Choice) {
    debug_assert_eq!(S, L + 1);

    let mut m = [0u64; L];
    let mut i = 0;
    while i < L {
        m[i] = msat[i];
        i += 1;
    }
    let m_inv = minv(m[0]);

    let mut f = *msat;
    let mut g = [0u64; S];
    let mut i = 0;
    while i < L {
        g[i] = a[i];
        i += 1;
    }

    let mut v = [0u64; L];
    let mut r = *one;
    // `eta` is `-delta`, and the algorithm starts at `delta = 1`
    let mut eta = -1i64;

    let batches = iterations.div_ceil(BATCH as usize);
    let mut b = 0;
    while b < batches {
        let (ne, t) = divsteps(eta, f[0], g[0]);

        let nf = lincomb_shift(t.ff, &f, t.fg, &g);
        let ng = lincomb_shift(t.gf, &f, t.gg, &g);
        let nv = update_acc::<L, S>(&m, m_inv, t.ff, &v, t.fg, &r);
        let nr = update_acc::<L, S>(&m, m_inv, t.gf, &v, t.gg, &r);

        eta = ne;
        f = nf;
        g = ng;
        v = nv;
        r = nr;
        b += 1;
    }

    // `g` is zero and `f` is ±1 by now, so the accumulator is `sign(f) / a`.
    // That convergence is what the step count is chosen to guarantee, so check
    // it in debug builds rather than trusting the bound silently.
    debug_assert!(
        converged(&f, &g),
        "safegcd did not converge in {} steps",
        iterations
    );

    (v, (f[S - 1] >> 63).ct_nonzero())
}

/// Whether the divsteps reached the terminal state: `g == 0` and `f == +/- 1`.
///
/// Only for `debug_assert!`
fn converged<const S: usize>(f: &[u64; S], g: &[u64; S]) -> bool {
    let mut i = 0;
    while i < S {
        if g[i] != 0 {
            return false;
        }
        i += 1;
    }
    // f is ±1 in two's complement: [1, 0, ..] or [!0, !0, ..]
    let neg = f[S - 1] >> 63 != 0;
    let (low, rest) = if neg { (!0u64, !0u64) } else { (1u64, 0u64) };
    if f[0] != low {
        return false;
    }
    let mut i = 1;
    while i < S {
        if f[i] != rest {
            return false;
        }
        i += 1;
    }
    true
}
