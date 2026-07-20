//! SM2 curve, defined over the prime field of order 2^256 - 2^224 - 2^96 + 2^64 - 1
//!
//! This is the recommended 256-bit prime-field curve of the Chinese SM2
//! standard (GB/T 32918 / GM/T 0003). It is a short Weierstrass curve with
//! `a = p - 3`, so it uses the same `a == -3` optimised point arithmetic as
//! the NIST/SEC2 `r1` curves.
use crate::curve::fiat::sm2_64::*;
use crate::curve::fiat::sm2_scalar_64::*;
use crate::curve::field::{Field, FieldSqrt, Sign};
use crate::curve::{
    affine, projective,
    weierstrass::{WeierstrassCurve, WeierstrassCurveAM3},
};
use crate::mp::ct::{Choice, CtEqual, CtOption, CtZero};
use crate::params::sm2::*;
use crate::{
    fiat_define_weierstrass_curve, fiat_define_weierstrass_curve_am3,
    fiat_define_weierstrass_points,
};
use crate::{fiat_field_montgomery_impl, fiat_field_sqrt_define};

const GM_LIMBS_SIZE: usize = 4;
const FE_LIMBS_SIZE: usize = 4;

fiat_field_montgomery_impl!(
    #[doc = "Element of the prime field Fp where p = 2^256 - 2^224 - 2^96 + 2^64 - 1"]
    FieldElement,
    256,
    P_LIMBS,
    FE_LIMBS_SIZE,
    fiat_sm2_non_montgomery_domain_field_element,
    fiat_sm2_nonzero,
    fiat_sm2_add,
    fiat_sm2_sub,
    fiat_sm2_mul,
    fiat_sm2_square,
    fiat_sm2_opp,
    fiat_sm2_to_bytes,
    fiat_sm2_from_bytes,
    fiat_sm2_montgomery_domain_field_element,
    fiat_sm2_to_montgomery,
    fiat_sm2_from_montgomery,
    fiat_sm2_selectznz,
    fiat_sm2_msat,
    fiat_sm2_divstep,
    fiat_sm2_divstep_precomp
);
fiat_field_sqrt_define!(FieldElement);

impl FieldElement {
    /// Get the multiplicative inverse
    ///
    /// Computes `self^(p-2)` via a repunit-based Fermat addition chain; the
    /// prime `p = 2^256 - 2^224 - 2^96 + 2^64 - 1` has long runs of set bits
    /// which the chain exploits. Cross-checked against the generic safegcd
    /// inversion by the test-suite.
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    pub fn inverse(&self) -> Self {
        assert!(!self.is_zero());
        let r1 = self.clone();
        let r2 = r1.square_rep(1) * &r1;
        let r4 = r2.square_rep(2) * &r2;
        let r8 = r4.square_rep(4) * &r4;
        let r16 = r8.square_rep(8) * &r8;
        let r32 = r16.square_rep(16) * &r16;
        let r64 = r32.square_rep(32) * &r32;
        let r128 = r64.square_rep(64) * &r64;
        let r3 = r2.square() * self;
        let r6 = r3.square_rep(3) * &r3;
        let r7 = r6.square() * self;
        let r14 = r7.square_rep(7) * &r7;
        let r15 = r14.square() * self;
        let r30 = r15.square_rep(15) * &r15;
        let r31 = r30.square() * self;
        let r62 = r31.square_rep(31) * &r31;

        let mut t1 = r31.clone();
        t1 = t1.square_rep(1);
        t1 = t1.square_rep(128) * &r128;
        t1 = t1.square_rep(32);
        t1 = t1.square_rep(62) * &r62;
        t1 = t1.square_rep(1);
        t1 = t1.square_rep(1) * &r1;
        t1
    }

    /// Compute the square root 'x' of the field element such that x*x = self
    ///
    /// Since `p = 3 (mod 4)`, the candidate root is `self^((p+1)/4)`; the
    /// returned `CtOption` is present only when it squares back to `self`.
    pub fn sqrt(&self) -> CtOption<Self> {
        let r1 = self.clone();
        let r2 = r1.square_rep(1) * &r1;
        let r4 = r2.square_rep(2) * &r2;
        let r8 = r4.square_rep(4) * &r4;
        let r16 = r8.square_rep(8) * &r8;
        let r32 = r16.square_rep(16) * &r16;
        let r64 = r32.square_rep(32) * &r32;
        let r128 = r64.square_rep(64) * &r64;
        let r3 = r2.square() * self;
        let r6 = r3.square_rep(3) * &r3;
        let r7 = r6.square() * self;
        let r14 = r7.square_rep(7) * &r7;
        let r15 = r14.square() * self;
        let r30 = r15.square_rep(15) * &r15;
        let r31 = r30.square() * self;

        let mut t1 = r31.clone();
        t1 = t1.square_rep(1);
        t1 = t1.square_rep(128) * &r128;
        t1 = t1.square_rep(31);
        t1 = t1.square_rep(1) * &r1;
        t1 = t1.square_rep(62);

        // t1 is the candidate root self^((p+1)/4); accept it only if it
        // squares back to self (i.e. self was actually a quadratic residue).
        let cand = t1;
        let cand_sq = &cand * &cand;
        CtOption::from((CtEqual::ct_eq(&cand_sq, self), cand))
    }
}

fiat_field_montgomery_impl!(
    #[doc = "Element of the prime field Fp for scalar where p is the order of the SM2 curve"]
    Scalar,
    256,
    ORDER_LIMBS,
    GM_LIMBS_SIZE,
    fiat_sm2_scalar_non_montgomery_domain_field_element,
    fiat_sm2_scalar_nonzero,
    fiat_sm2_scalar_add,
    fiat_sm2_scalar_sub,
    fiat_sm2_scalar_mul,
    fiat_sm2_scalar_square,
    fiat_sm2_scalar_opp,
    fiat_sm2_scalar_to_bytes,
    fiat_sm2_scalar_from_bytes,
    fiat_sm2_scalar_montgomery_domain_field_element,
    fiat_sm2_scalar_to_montgomery,
    fiat_sm2_scalar_from_montgomery,
    fiat_sm2_scalar_selectznz,
    fiat_sm2_scalar_msat,
    fiat_sm2_scalar_divstep,
    fiat_sm2_scalar_divstep_precomp
);

impl Scalar {
    /// Get the multiplicative inverse
    ///
    /// Computes `self^(n-2)` via a Fermat addition chain, where `n` is the
    /// group order. Cross-checked against the generic safegcd inversion by
    /// the test-suite.
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    /// TODO this will change to being a method of NonZeroScalar
    pub fn inverse(&self) -> Self {
        assert!(!self.is_zero());
        let r1 = self.clone();
        let r2 = r1.square_rep(1) * &r1;
        let r3 = r2.square() * self;
        let r6 = r3.square_rep(3) * &r3;
        let r12 = r6.square_rep(6) * &r6;
        let r24 = r12.square_rep(12) * &r12;
        let r48 = r24.square_rep(24) * &r24;
        let r96 = r48.square_rep(48) * &r48;
        let r7 = r6.square() * self;
        let r14 = r7.square_rep(7) * &r7;
        let r15 = r14.square() * self;
        let r30 = r15.square_rep(15) * &r15;
        let r31 = r30.square() * self;
        let r4 = r2.square_rep(2) * &r2;
        let r5 = r4.square() * self;

        let mut t = r31.clone();
        t = t.square_rep(1);
        t = t.square_rep(96) * &r96;
        t = t.square_rep(1);
        t = t.square_rep(3) * &r3;
        t = t.square_rep(2);
        t = t.square_rep(1) * &r1;
        t = t.square_rep(7);
        t = t.square_rep(4) * &r4;
        t = t.square_rep(1);
        t = t.square_rep(5) * &r5;
        t = t.square_rep(1);
        t = t.square_rep(2) * &r2;
        t = t.square_rep(1);
        t = t.square_rep(1) * &r1;
        t = t.square_rep(1);
        t = t.square_rep(2) * &r2;
        t = t.square_rep(2);
        t = t.square_rep(1) * &r1;
        t = t.square_rep(4);
        t = t.square_rep(3) * &r3;
        t = t.square_rep(3);
        t = t.square_rep(2) * &r2;
        t = t.square_rep(6);
        t = t.square_rep(1) * &r1;
        t = t.square_rep(1);
        t = t.square_rep(1) * &r1;
        t = t.square_rep(2);
        t = t.square_rep(1) * &r1;
        t = t.square_rep(1);
        t = t.square_rep(1) * &r1;
        t = t.square_rep(1);
        t = t.square_rep(2) * &r2;
        t = t.square_rep(1);
        t = t.square_rep(1) * &r1;
        t = t.square_rep(1);
        t = t.square_rep(1) * &r1;
        t = t.square_rep(2);
        t = t.square_rep(3) * &r3;
        t = t.square_rep(1);
        t = t.square_rep(3) * &r3;
        t = t.square_rep(1);
        t = t.square_rep(6) * &r6;
        t = t.square_rep(1);
        t = t.square_rep(1) * &r1;
        t = t.square_rep(6);
        t = t.square_rep(1) * &r1;
        t = t.square_rep(2);
        t = t.square_rep(1) * &r1;
        t = t.square_rep(2);
        t = t.square_rep(3) * &r3;
        t = t.square_rep(2);
        t = t.square_rep(3) * &r3;
        t = t.square_rep(1);
        t = t.square_rep(1) * &r1;
        t = t.square_rep(1);
        t = t.square_rep(1) * &r1;
        t = t.square_rep(1);
        t = t.square_rep(1) * &r1;
        t = t.square_rep(1);
        t = t.square_rep(1) * &r1;
        t = t.square_rep(5);
        t = t.square_rep(1) * &r1;
        t = t.square_rep(2);
        t = t.square_rep(1) * &r1;
        t = t.square_rep(4);
        t = t.square_rep(1) * &r1;
        t
    }
}

fiat_define_weierstrass_curve!(FieldElement);
fiat_define_weierstrass_curve_am3!(FieldElement);
fiat_define_weierstrass_points!(FieldElement);

#[cfg(test)]
mod tests {
    mod fe {
        use super::super::FieldElement;
        use crate::{fiat_field_sqrt_unittest, fiat_field_unittest};

        fiat_field_unittest!(FieldElement);
        crate::fiat_field_safegcd_unittest!(FieldElement);
        fiat_field_sqrt_unittest!(FieldElement);
    }
    mod gm {
        use super::super::Scalar;
        use crate::fiat_field_unittest;
        fiat_field_unittest!(Scalar);
        crate::fiat_field_safegcd_unittest!(Scalar);
    }

    mod endian {
        use super::super::{FieldElement, Scalar};

        #[test]
        fn default_endianness_is_be() {
            // SM2 uses big-endian, like the SEC2 curves
            let x = FieldElement::from_u64(0x1234_5678_9abc_def0);
            assert_eq!(x.to_bytes(), x.to_bytes_be());
            assert_ne!(x.to_bytes(), x.to_bytes_le());
            assert_eq!(FieldElement::from_bytes(&x.to_bytes_be()).unwrap(), x);

            let s = Scalar::from_u64(0x0fed_cba9_8765_4321);
            assert_eq!(s.to_bytes(), s.to_bytes_be());
            assert_eq!(Scalar::from_bytes(&s.to_bytes_be()).unwrap(), s);
        }
    }

    mod point {
        crate::fiat_curve_point_unittest!();
    }

    mod group {
        use super::super::{Point, PointAffine, Scalar};
        use crate::curve::group::CurveGroup;

        /// The generator must lie on the curve.
        #[test]
        fn generator_on_curve() {
            let (gx, gy) = (PointAffine::GENERATOR).to_coordinate();
            assert!(PointAffine::from_coordinate(gx, gy).is_some());
        }

        /// The generator has order `n` (the group order), so `n * G = O`.
        ///
        /// `n` itself is not a representable scalar (scalars live in `[0, n)`),
        /// but `-1 (mod n) = n - 1` is. Multiplying the generator by it must
        /// give `-G`, since `(n - 1) * G = n*G - G = -G` exactly when `n*G = O`.
        #[test]
        fn generator_has_group_order() {
            let n_minus_1 = -Scalar::one();
            let r = &Point::GENERATOR * &n_minus_1;
            assert_eq!(r, -&Point::GENERATOR);
            assert_eq!(&r + &Point::GENERATOR, Point::INFINITY);
        }

        /// The fixed-base comb (`mul_base`) agrees with the general scalar
        /// multiplication of the generator.
        #[test]
        fn mul_base_matches_generic() {
            for k in [1u64, 2, 3, 7, 0x1234_5678, 0xdead_beef_cafe] {
                let s = Scalar::from_u64(k);
                let base = Point::mul_base(&s);
                let generic = &Point::GENERATOR * &s;
                assert_eq!(base, generic, "mul_base != G * k for k={}", k);
            }
        }

        /// `mul_vartime` agrees with the constant-time scalar multiplication.
        #[test]
        fn mul_vartime_matches_ct() {
            let two = Scalar::from_u64(2);
            let g2 = Point::GENERATOR.double();
            assert_eq!(Point::GENERATOR.mul_vartime(&two), g2);
            assert_eq!(&Point::GENERATOR * &two, g2);
        }

        /// Known-answer tests for `k * G`, computed with an independent
        /// reference implementation of the SM2 curve. This pins
        /// the whole stack (field arithmetic, point add/double, and the
        /// embedded generator) to the real SM2 curve.
        #[test]
        fn scalar_mul_kats() {
            use super::super::FieldElement;
            let kats: [(u64, [u8; 32], [u8; 32]); 3] = [
                (
                    2,
                    [
                        0x56, 0xce, 0xfd, 0x60, 0xd7, 0xc8, 0x7c, 0x00, 0x0d, 0x58, 0xef, 0x57,
                        0xfa, 0x73, 0xba, 0x4d, 0x9c, 0x0d, 0xfa, 0x08, 0xc0, 0x8a, 0x73, 0x31,
                        0x49, 0x5c, 0x2e, 0x1d, 0xa3, 0xf2, 0xbd, 0x52,
                    ],
                    [
                        0x31, 0xb7, 0xe7, 0xe6, 0xcc, 0x81, 0x89, 0xf6, 0x68, 0x53, 0x5c, 0xe0,
                        0xf8, 0xea, 0xf1, 0xbd, 0x6d, 0xe8, 0x4c, 0x18, 0x2f, 0x6c, 0x8e, 0x71,
                        0x6f, 0x78, 0x0d, 0x3a, 0x97, 0x0a, 0x23, 0xc3,
                    ],
                ),
                (
                    0x1_2345_6789,
                    [
                        0x96, 0xfc, 0xd2, 0x5f, 0x7e, 0x9e, 0xcb, 0xfa, 0xf7, 0xae, 0x7e, 0x39,
                        0x9d, 0x11, 0x17, 0xfc, 0x38, 0xc6, 0x73, 0x6e, 0x8f, 0xfd, 0x64, 0x0b,
                        0xa4, 0xf5, 0xcd, 0x37, 0x8a, 0x34, 0xff, 0xf4,
                    ],
                    [
                        0x99, 0x5a, 0xee, 0x11, 0x91, 0x6e, 0x80, 0x65, 0xb3, 0x8d, 0xa9, 0xec,
                        0xb0, 0xfd, 0x8b, 0xa4, 0x52, 0x80, 0xb3, 0x43, 0xd1, 0x39, 0x8a, 0xe6,
                        0x75, 0xeb, 0xcf, 0x09, 0xa5, 0xa0, 0x5d, 0x99,
                    ],
                ),
                (
                    0xdead_beef_cafe_1234,
                    [
                        0x2e, 0xaa, 0xe9, 0x2e, 0x05, 0x0d, 0x06, 0x64, 0x15, 0x46, 0x75, 0x8c,
                        0xf8, 0xbc, 0x92, 0xb5, 0x4e, 0xd4, 0x33, 0xe6, 0x74, 0x05, 0xd0, 0xc5,
                        0x46, 0xdf, 0xd2, 0x7b, 0x56, 0xdf, 0xb8, 0x7d,
                    ],
                    [
                        0x6b, 0x28, 0xd9, 0x29, 0x17, 0x10, 0x31, 0x6e, 0xf9, 0x7a, 0xcd, 0x3a,
                        0x74, 0x1e, 0x34, 0x07, 0x3f, 0x3a, 0x4d, 0x26, 0xef, 0x0f, 0x28, 0xbe,
                        0x84, 0x7b, 0x01, 0x41, 0x93, 0xbc, 0x00, 0x03,
                    ],
                ),
            ];
            for (k, ex, ey) in kats {
                let s = Scalar::from_u64(k);
                let ex = FieldElement::from_bytes(&ex).unwrap();
                let ey = FieldElement::from_bytes(&ey).unwrap();
                for p in [&Point::GENERATOR * &s, Point::mul_base(&s)] {
                    let aff = p.to_affine().expect("k*G is not infinity");
                    let (x, y) = aff.to_coordinate();
                    assert_eq!(x, &ex, "x mismatch for k={:#x}", k);
                    assert_eq!(y, &ey, "y mismatch for k={:#x}", k);
                }
            }
        }
    }
}
