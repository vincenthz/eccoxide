//! BLS12-381 scalar field `Scalar` (a.k.a. `Fr`)
//!
//! `Scalar` is the prime field of order `r`, the 255-bit order of the
//! prime-order subgroups G1 (and G2). It is the field of scalars acting on the
//! curve group, so scalar multiplication of a point takes a `Scalar`.
//!
//! As for [`super::fp`], the arithmetic comes from the fiat-crypto Montgomery
//! backend and this module only wires it into the generic [`Field`]
//! abstraction.

use crate::curve::fiat::bls12_381_scalar_64::*;
use crate::curve::field::{Field, Sign};
use crate::fiat_field_montgomery_impl;
use crate::mp::ct::{Choice, CtEqual, CtOption, CtZero};
use crate::params::bls12_381::ORDER_LIMBS;

/// Number of 64-bit limbs of a scalar field element (255 bits -> 4 limbs).
const GM_LIMBS_SIZE: usize = 4;

fiat_field_montgomery_impl!(
    #[doc = "Element of the BLS12-381 scalar prime field Fr where r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001"]
    Scalar,
    255,
    ORDER_LIMBS,
    GM_LIMBS_SIZE,
    fiat_bls12_381_scalar_non_montgomery_domain_field_element,
    fiat_bls12_381_scalar_nonzero,
    fiat_bls12_381_scalar_add,
    fiat_bls12_381_scalar_sub,
    fiat_bls12_381_scalar_mul,
    fiat_bls12_381_scalar_square,
    fiat_bls12_381_scalar_opp,
    fiat_bls12_381_scalar_to_bytes,
    fiat_bls12_381_scalar_from_bytes,
    fiat_bls12_381_scalar_montgomery_domain_field_element,
    fiat_bls12_381_scalar_to_montgomery,
    fiat_bls12_381_scalar_from_montgomery,
    fiat_bls12_381_scalar_selectznz,
    fiat_bls12_381_scalar_msat,
    fiat_bls12_381_scalar_divstep,
    fiat_bls12_381_scalar_divstep_precomp
);

// Tonelli-Shanks square-root constants for `r`.
//
// `r - 1 = t * 2^S` with `t` odd and `S = 32` (the 2-adicity of the field).
// `ROOT_OF_UNITY` is a generator of the order-`2^S` subgroup, taken as `z^t`
// for the least quadratic non-residue `z = 5`. All big-endian.
const SQRT_S: u32 = 32;
/// `t = (r - 1) / 2^S`, the odd part of `r - 1` (big-endian).
const SQRT_T_BYTES: [u8; 32] = [
    0x00, 0x00, 0x00, 0x00, 0x73, 0xed, 0xa7, 0x53, 0x29, 0x9d, 0x7d, 0x48, 0x33, 0x39, 0xd8, 0x08,
    0x09, 0xa1, 0xd8, 0x05, 0x53, 0xbd, 0xa4, 0x02, 0xff, 0xfe, 0x5b, 0xfe, 0xff, 0xff, 0xff, 0xff,
];
/// `(t + 1) / 2`, the exponent of the initial candidate root (big-endian).
const SQRT_T_PLUS1_DIV2_BYTES: [u8; 32] = [
    0x00, 0x00, 0x00, 0x00, 0x39, 0xf6, 0xd3, 0xa9, 0x94, 0xce, 0xbe, 0xa4, 0x19, 0x9c, 0xec, 0x04,
    0x04, 0xd0, 0xec, 0x02, 0xa9, 0xde, 0xd2, 0x01, 0x7f, 0xff, 0x2d, 0xff, 0x80, 0x00, 0x00, 0x00,
];
/// `ROOT_OF_UNITY = z^t` (`z = 5`), a primitive `2^S`-th root of unity (big-endian).
const SQRT_ROOT_OF_UNITY_BYTES: [u8; 32] = [
    0x02, 0x12, 0xd7, 0x9e, 0x5b, 0x41, 0x6b, 0x6f, 0x0f, 0xd5, 0x6d, 0xc8, 0xd1, 0x68, 0xd6, 0xc0,
    0xc4, 0x02, 0x4f, 0xf2, 0x70, 0xb3, 0xe0, 0x94, 0x1b, 0x78, 0x8f, 0x50, 0x0b, 0x91, 0x2f, 0x1f,
];

impl Scalar {
    /// Get the multiplicative inverse.
    ///
    /// This currently delegates to the representation-agnostic Bernstein-Yang
    /// "safegcd" inversion ([`Self::inverse_safegcd`]); see the note on
    /// [`super::fp::Fp::inverse`].
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a
    /// panic.
    pub fn inverse(&self) -> Self {
        self.inverse_safegcd()
    }

    /// Compute a square root `x` such that `x*x = self`, or `None` when `self`
    /// is not a quadratic residue.
    ///
    /// Unlike [`super::fp::Fp::sqrt`], `r ≡ 1 (mod 2^32)` (2-adicity 32), so
    /// there is no simple exponentiation formula; this uses the Tonelli-Shanks
    /// algorithm. The returned root is the one produced by the algorithm; the
    /// two roots are `x` and `-x`.
    ///
    /// This is used for (de)compressing points of the Jubjub curve (whose base
    /// field is this scalar field). Note that the running time depends on the
    /// input (the number of Tonelli-Shanks iterations is value-dependent), so
    /// it is **not** constant-time; only use it on public values.
    pub fn sqrt(&self) -> CtOption<Self> {
        // sqrt(0) = 0
        if self.is_zero() {
            return CtOption::from((Choice::TRUE, Self::zero()));
        }

        let one = Self::one();
        // candidate root and the "tracking" element t = self^t_odd
        let mut r = self.power(&SQRT_T_PLUS1_DIV2_BYTES);
        let mut t = self.power(&SQRT_T_BYTES);
        let mut c = Self::from_bytes_unchecked_be(&SQRT_ROOT_OF_UNITY_BYTES);
        let mut m = SQRT_S;

        while t != one {
            // least i in 1..m with t^(2^i) == 1
            let mut i = 0u32;
            let mut t2i = t.clone();
            while t2i != one {
                t2i = t2i.square();
                i += 1;
                if i == m {
                    // t has full order 2^m: self is a non-residue
                    return CtOption::from((Choice::FALSE, Self::zero()));
                }
            }
            // b = c^(2^(m - i - 1)); here i < m so the exponent is non-negative
            let mut b = c.clone();
            for _ in 0..(m - i - 1) {
                b = b.square();
            }
            r = &r * &b;
            c = b.square();
            t = &t * &c;
            m = i;
        }

        // defensive re-check (always true for a residue once we reach here)
        let present = CtEqual::ct_eq(&r.square(), self);
        CtOption::from((present, r))
    }
}

#[cfg(test)]
mod tests {
    use super::Scalar;
    use crate::fiat_field_unittest;

    fiat_field_unittest!(Scalar);
    crate::fiat_field_safegcd_unittest!(Scalar);

    #[test]
    fn sqrt_roundtrip() {
        // for every square s = x^2, sqrt(s)^2 == s (the root is x or -x)
        for k in [1u64, 2, 3, 4, 5, 7, 100, 123456789, 0x0123_4567_89ab_cdef] {
            let x = Scalar::from_u64(k);
            let s = x.square();
            let root = s.sqrt().into_option().expect("square has a root");
            assert_eq!(root.square(), s, "k={}", k);
            assert!(root == x || root == -&x, "k={}", k);
        }
        // sqrt(0) == 0
        assert_eq!(Scalar::zero().sqrt().into_option(), Some(Scalar::zero()));
    }

    #[test]
    fn sqrt_non_residue() {
        // 5 is the least quadratic non-residue mod r
        assert!(Scalar::from_u64(5).sqrt().into_option().is_none());
    }
}
