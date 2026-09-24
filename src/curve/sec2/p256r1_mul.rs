//! Montgomery multiplication specialised for the p256 base field.
//!
//! Hand-written replacement for the fiat-crypto generated `fiat_p256_mul` and
//! `fiat_p256_square`, exploiting the shape of
//!
//! ```text
//! p = 2^256 - 2^224 + 2^192 + 2^96 - 1
//!   = 0xffffffff00000001_000000000000000_000000000fffffff_fffffffffffffffff
//!   = [0xffffffffffffffff, 0x00000000ffffffff, 0, 0xffffffff00000001]
//! ```
//!
//! The schoolbook multiplication part is unchanged, only the reduction is
//! optimised for the special form of `p`
//!
//! # The reduction
//!
//! `p = -1 (mod 2^64)`, so `-p^-1 = 1 (mod 2^64)` and the Montgomery multiplier
//! of a reduction round is just the low limb of the accumulator, `m = t[0]`,
//! with no multiplication to derive it. Adding `m * p` to the accumulator and
//! shifting one limb down is then
//!
//! ```text
//! (t + m*p) / 2^64 = (t - m)/2^64 + m*2^32 + m*(2^64 - 2^32 + 1)*2^128
//! ```
//!
//! using `m*p = m*2^256 - m*2^224 + m*2^192 + m*2^96 - m`. The three parts are
//!
//! * `(t - m)/2^64`: an exact shift, since `m` was chosen to cancel `t[0]`;
//! * `m * 2^32`, contributing `m << 32` to limb 0 and `m >> 32` to limb 1;
//! * `m * 0xffffffff00000001` at limb offset 2, and since that constant is
//!   `2^64 - 2^32 + 1` the 128-bit product is `(m - (m>>32))*2^64 + (m - (m<<32))`
//!   with the borrow of the low half folded into the high half.
//!
//! A reduction round is:
//!
//! * two shifts
//! * a two-limb subtract chain
//! * a four-limb add chain
//!
//! no word multiplication at all, against the 12 the generic word-by-word Montgomery code issues
//! (three per round, one per non-trivial limb of `p`).
//!
//! Both bounds need `a` and `b` reduced

#[cfg(not(feature = "p256r1-optimised"))]
use crate::curve::fiat::p256_64::fiat_p256_montgomery_domain_field_element;

#[cfg(feature = "p256r1-optimised")]
mod optimised {
    use crate::curve::fiat::p256_64::fiat_p256_montgomery_domain_field_element;

    /// `p`, little-endian limbs.
    const P: [u64; 4] = [
        0xffff_ffff_ffff_ffff,
        0x0000_0000_ffff_ffff,
        0x0000_0000_0000_0000,
        0xffff_ffff_0000_0001,
    ];

    /// Widening multiply, returning `(low, high)`.
    #[inline(always)]
    const fn wide(a: u64, b: u64) -> (u64, u64) {
        let t = (a as u128) * (b as u128);
        (t as u64, (t >> 64) as u64)
    }

    /// `a + b + carry`, returning the sum and the carry out.
    #[inline(always)]
    const fn adc(a: u64, b: u64, carry: bool) -> (u64, bool) {
        let (x, c1) = a.overflowing_add(b);
        let (y, c2) = x.overflowing_add(carry as u64);
        // at most one of c1, c2 can be set
        (y, c1 | c2)
    }

    /// `a - b - borrow`, returning the difference and the borrow out.
    #[inline(always)]
    const fn sbb(a: u64, b: u64, borrow: bool) -> (u64, bool) {
        let (x, b1) = a.overflowing_sub(b);
        let (y, b2) = x.overflowing_sub(borrow as u64);
        // at most one of b1, b2 can be set
        (y, b1 | b2)
    }

    /// `a * b * 2^-256 mod p`, for `a` and `b` already reduced.
    #[inline(always)]
    const fn mont_mul(a: &[u64; 4], b: &[u64; 4]) -> [u64; 4] {
        let (a0, a1, a2, a3) = (a[0], a[1], a[2], a[3]);
        let (mut r0, mut r1, mut r2, mut r3, mut r4) = (0u64, 0u64, 0u64, 0u64, 0u64);

        let mut i = 0;
        while i < 4 {
            let bi = b[i];

            let (l0, h0) = wide(a0, bi);
            let (l1, h1) = wide(a1, bi);
            let (l2, h2) = wide(a2, bi);
            let (l3, h3) = wide(a3, bi);

            // r += a * b[i], as two passes so that each is a single carry chain:
            // first the low halves of the four products, then the high halves
            // shifted up by one limb. Neither carries out of the fifth limb, the
            // sum being below 2^320 (see the bounds above).
            let (x0, c) = adc(r0, l0, false);
            let (x1, c) = adc(r1, l1, c);
            let (x2, c) = adc(r2, l2, c);
            let (x3, c) = adc(r3, l3, c);
            let (x4, _) = adc(r4, 0, c);

            let (y1, c) = adc(x1, h0, false);
            let (y2, c) = adc(x2, h1, c);
            let (y3, c) = adc(x3, h2, c);
            let (y4, _) = adc(x4, h3, c);

            // r = (r + m*p) / 2^64 with m = r[0]; see the module documentation.
            // The second pass starts at limb 1, so limb 0 is final after the first.
            let m = x0;
            let lo32 = m << 32;
            let hi32 = m >> 32;
            // d = m * 0xffffffff00000001 = m*2^64 - m*2^32 + m, which is
            // (m - (m>>32))*2^64 + (m - (m<<32)) and never borrows out, since
            // m >= (m>>32) + 1 whenever the low half does borrow.
            let (d0, brw) = sbb(m, lo32, false);
            let (d1, _) = sbb(m, hi32, brw);

            let (z0, c) = adc(y1, lo32, false);
            let (z1, c) = adc(y2, hi32, c);
            let (z2, c) = adc(y3, d0, c);
            let (z3, c) = adc(y4, d1, c);

            r0 = z0;
            r1 = z1;
            r2 = z2;
            r3 = z3;
            r4 = c as u64;

            i += 1;
        }

        // r < 2p: subtract p once if that does not borrow.
        let (s0, brw) = sbb(r0, P[0], false);
        let (s1, brw) = sbb(r1, P[1], brw);
        let (s2, brw) = sbb(r2, P[2], brw);
        let (s3, brw) = sbb(r3, P[3], brw);
        let (_, brw) = sbb(r4, 0, brw);

        let keep = 0u64.wrapping_sub(brw as u64);
        [
            (r0 & keep) | (s0 & !keep),
            (r1 & keep) | (s1 & !keep),
            (r2 & keep) | (s2 & !keep),
            (r3 & keep) | (s3 & !keep),
        ]
    }

    /// Montgomery multiplication in the p256 base field.
    ///
    /// Drop-in replacement for `fiat_p256_mul`: `out1 = arg1 * arg2 * 2^-256 mod p`,
    /// with both arguments and the result in the Montgomery domain and reduced.
    pub const fn p256_mul(
        out1: &mut fiat_p256_montgomery_domain_field_element,
        arg1: &fiat_p256_montgomery_domain_field_element,
        arg2: &fiat_p256_montgomery_domain_field_element,
    ) {
        out1.0 = mont_mul(&arg1.0, &arg2.0);
    }

    /// Montgomery squaring in the p256 base field.
    ///
    /// Drop-in replacement for `fiat_p256_square`. Like the fiat routine it reuses
    /// the multiplication schedule rather than a symmetry-exploiting one, but since
    /// [`mont_mul`] is inlined here the backend does get to share the `a[i]*a[j]`
    /// and `a[j]*a[i]` products, which brings the 16 word products down to 10.
    pub const fn p256_square(
        out1: &mut fiat_p256_montgomery_domain_field_element,
        arg1: &fiat_p256_montgomery_domain_field_element,
    ) {
        out1.0 = mont_mul(&arg1.0, &arg1.0);
    }
}

#[cfg(feature = "p256r1-optimised")]
pub use optimised::{p256_mul, p256_square};

#[cfg(not(feature = "p256r1-optimised"))]
#[inline]
pub const fn p256_mul(
    out1: &mut fiat_p256_montgomery_domain_field_element,
    arg1: &fiat_p256_montgomery_domain_field_element,
    arg2: &fiat_p256_montgomery_domain_field_element,
) {
    use super::super::fiat::p256_64::fiat_p256_mul;
    fiat_p256_mul(out1, arg1, arg2)
}

#[cfg(not(feature = "p256r1-optimised"))]
#[inline]
pub const fn p256_square(
    out1: &mut fiat_p256_montgomery_domain_field_element,
    arg1: &fiat_p256_montgomery_domain_field_element,
) {
    use super::super::fiat::p256_64::fiat_p256_square;
    fiat_p256_square(out1, arg1)
}
