//! BLS12-381 G1
//!
//! G1 is the prime-order-`r` subgroup of the short Weierstrass curve
//! `E(Fp): y^2 = x^3 + 4` (so `a = 0`, `b = 4`) defined over the base field
//! [`Fp`](super::fp::Fp). Its scalar field is [`Scalar`](super::scalar::Scalar).
//!
//! This wires the generic curve machinery (affine / projective points, the
//! complete addition formulas and the fixed-base comb) into the [`CurveGroup`]
//! abstraction.
//!
//! [`CurveGroup`]: crate::curve::group::CurveGroup

use crate::curve::bls12_381::fp::Fp as FieldElement;
use crate::curve::bls12_381::scalar::Scalar;
use crate::curve::field::{Field, Sign};
use crate::curve::{
    affine, projective,
    weierstrass::{WeierstrassCurve, WeierstrassCurveA0},
};
use crate::params::bls12_381::*;
use crate::{
    bls12_381_define_point_serialization, fiat_define_weierstrass_curve,
    fiat_define_weierstrass_curve_a0, fiat_define_weierstrass_points,
};

fiat_define_weierstrass_curve!(FieldElement);
fiat_define_weierstrass_curve_a0!(FieldElement);
fiat_define_weierstrass_points!(FieldElement);
bls12_381_define_point_serialization!(FieldElement);

#[cfg(test)]
mod tests {
    use super::{FieldElement, Point, PointAffine, Scalar};
    use crate::curve::group::CurveGroup;

    fn s(v: u64) -> Scalar {
        Scalar::from_u64(v)
    }

    #[test]
    fn generator_on_curve() {
        // the generator coordinates must satisfy the curve equation
        let (gx, gy) = PointAffine::GENERATOR.to_coordinate();
        assert!(PointAffine::from_coordinate(gx, gy).is_some());
    }

    #[test]
    fn identity_laws() {
        let g = Point::GENERATOR;
        let id = Point::INFINITY;
        assert_eq!(&g + &id, g);
        assert_eq!(&id + &g, g);
        assert_eq!(&g - &g, id);
        assert_eq!(&g + (-&g), id);
    }

    #[test]
    fn double_is_add() {
        let g = Point::GENERATOR;
        assert_eq!(g.double(), &g + &g);
    }

    #[test]
    fn scalar_mul_matches_addition() {
        let g = Point::GENERATOR;
        let g2 = g.double();
        let g3 = &g2 + &g;
        let three = s(3);
        assert_eq!(&g * &three, g3);
        assert_eq!(Point::mul_base(&three), g3);
        assert_eq!(g.mul_vartime(&three), g3);
    }

    #[test]
    fn mul_base_matches_variable_base() {
        // the fixed-base comb table must agree with variable-base multiplication
        let g = Point::GENERATOR;
        for v in [1u64, 2, 7, 255, 256, 0x1_0001, 0x1234_5678_9abc_def0] {
            let k = s(v);
            assert_eq!(Point::mul_base(&k), &g * &k, "mul_base != g*k for {}", v);
            assert_eq!(g.mul_vartime(&k), &g * &k, "mul_vartime != g*k for {}", v);
        }
    }

    #[test]
    fn scalar_mul_is_linear() {
        let g = Point::GENERATOR;
        let a = Scalar::from(0x1234_5678_9abc_def0u64);
        let b = Scalar::from(0x0fed_cba9_8765_4321u64);
        let ab = &a + &b;
        assert_eq!(
            Point::mul_base(&ab),
            Point::mul_base(&a) + Point::mul_base(&b)
        );
        assert_eq!(&g * &ab, &g * &a + &g * &b);
    }

    #[test]
    fn compress_decompress_roundtrip() {
        let p = (Point::GENERATOR * &s(12345)).to_affine().unwrap();
        let (x, sign) = p.compress();
        let x = x.clone();
        let recovered = PointAffine::decompress(&x, sign).into_option().unwrap();
        assert_eq!(recovered, p);
    }

    #[test]
    fn to_affine_from_affine_roundtrip() {
        let p = Point::GENERATOR * &s(98765);
        let affine = p.to_affine().unwrap();
        assert_eq!(Point::from_affine(&affine), p);
    }

    #[test]
    fn generator_field_element_bytes_roundtrip() {
        let (gx, _) = PointAffine::GENERATOR.to_coordinate();
        let bytes = gx.to_bytes();
        assert_eq!(&FieldElement::from_bytes(&bytes).unwrap(), gx);
    }

    crate::bls12_381_point_serialization_unittest!();

    /// Known-answer vectors for the standard G1 encodings: `k * G` for small
    /// `k`, byte for byte.
    ///
    /// `k = 1` is the published encoding of the generator: its x-coordinate
    /// `0x17f1d3a7...` with the compression flag set, and the sort flag clear
    /// since its y is below `(p - 1)/2`. `k = 2` and `k = 0xdeadbeef` have the
    /// sort flag set instead, so the vectors pin the flag both ways. The
    /// multiples were computed independently.
    #[test]
    fn serialization_kat() {
        const COMPRESSED: &[(u64, [u8; 48])] = &[
            (
                1,
                [
                    0x97, 0xf1, 0xd3, 0xa7, 0x31, 0x97, 0xd7, 0x94, 0x26, 0x95, 0x63, 0x8c, 0x4f,
                    0xa9, 0xac, 0x0f, 0xc3, 0x68, 0x8c, 0x4f, 0x97, 0x74, 0xb9, 0x05, 0xa1, 0x4e,
                    0x3a, 0x3f, 0x17, 0x1b, 0xac, 0x58, 0x6c, 0x55, 0xe8, 0x3f, 0xf9, 0x7a, 0x1a,
                    0xef, 0xfb, 0x3a, 0xf0, 0x0a, 0xdb, 0x22, 0xc6, 0xbb,
                ],
            ),
            (
                2,
                [
                    0xa5, 0x72, 0xcb, 0xea, 0x90, 0x4d, 0x67, 0x46, 0x88, 0x08, 0xc8, 0xeb, 0x50,
                    0xa9, 0x45, 0x0c, 0x97, 0x21, 0xdb, 0x30, 0x91, 0x28, 0x01, 0x25, 0x43, 0x90,
                    0x2d, 0x0a, 0xc3, 0x58, 0xa6, 0x2a, 0xe2, 0x8f, 0x75, 0xbb, 0x8f, 0x1c, 0x7c,
                    0x42, 0xc3, 0x9a, 0x8c, 0x55, 0x29, 0xbf, 0x0f, 0x4e,
                ],
            ),
            (
                3,
                [
                    0x89, 0xec, 0xe3, 0x08, 0xf9, 0xd1, 0xf0, 0x13, 0x17, 0x65, 0x21, 0x2d, 0xec,
                    0xa9, 0x96, 0x97, 0xb1, 0x12, 0xd6, 0x1f, 0x9b, 0xe9, 0xa5, 0xf1, 0xf3, 0x78,
                    0x0a, 0x51, 0x33, 0x5b, 0x3f, 0xf9, 0x81, 0x74, 0x7a, 0x0b, 0x2c, 0xa2, 0x17,
                    0x9b, 0x96, 0xd2, 0xc0, 0xc9, 0x02, 0x4e, 0x52, 0x24,
                ],
            ),
            (
                12345,
                [
                    0x85, 0x30, 0xc1, 0xbd, 0xc4, 0xcd, 0x6b, 0x14, 0x08, 0xbe, 0x09, 0x33, 0xc4,
                    0xa4, 0x1a, 0xc3, 0x51, 0x33, 0x50, 0xee, 0xf3, 0x68, 0x50, 0xb8, 0x04, 0x70,
                    0x8e, 0x1f, 0x33, 0x89, 0x32, 0xce, 0x01, 0xb6, 0x55, 0xa1, 0x63, 0x34, 0x4a,
                    0x45, 0x00, 0xb2, 0x81, 0xc8, 0x75, 0x0c, 0x46, 0x1f,
                ],
            ),
            (
                0xdead_beef,
                [
                    0xac, 0xcc, 0xb5, 0xba, 0xb2, 0x94, 0x4a, 0x1b, 0xdc, 0x72, 0x1c, 0x97, 0xf3,
                    0xaf, 0xfa, 0x03, 0x5d, 0x50, 0x7c, 0x78, 0xfe, 0x44, 0x2a, 0x92, 0x84, 0x98,
                    0x2b, 0xd4, 0xc2, 0x76, 0x17, 0xb3, 0x3f, 0x1d, 0x46, 0xe8, 0x19, 0x1a, 0x1e,
                    0xda, 0x03, 0xd7, 0x3c, 0x35, 0x77, 0x52, 0xd2, 0x19,
                ],
            ),
        ];
        const UNCOMPRESSED: &[(u64, [u8; 96])] = &[
            (
                1,
                [
                    0x17, 0xf1, 0xd3, 0xa7, 0x31, 0x97, 0xd7, 0x94, 0x26, 0x95, 0x63, 0x8c, 0x4f,
                    0xa9, 0xac, 0x0f, 0xc3, 0x68, 0x8c, 0x4f, 0x97, 0x74, 0xb9, 0x05, 0xa1, 0x4e,
                    0x3a, 0x3f, 0x17, 0x1b, 0xac, 0x58, 0x6c, 0x55, 0xe8, 0x3f, 0xf9, 0x7a, 0x1a,
                    0xef, 0xfb, 0x3a, 0xf0, 0x0a, 0xdb, 0x22, 0xc6, 0xbb, 0x08, 0xb3, 0xf4, 0x81,
                    0xe3, 0xaa, 0xa0, 0xf1, 0xa0, 0x9e, 0x30, 0xed, 0x74, 0x1d, 0x8a, 0xe4, 0xfc,
                    0xf5, 0xe0, 0x95, 0xd5, 0xd0, 0x0a, 0xf6, 0x00, 0xdb, 0x18, 0xcb, 0x2c, 0x04,
                    0xb3, 0xed, 0xd0, 0x3c, 0xc7, 0x44, 0xa2, 0x88, 0x8a, 0xe4, 0x0c, 0xaa, 0x23,
                    0x29, 0x46, 0xc5, 0xe7, 0xe1,
                ],
            ),
            (
                12345,
                [
                    0x05, 0x30, 0xc1, 0xbd, 0xc4, 0xcd, 0x6b, 0x14, 0x08, 0xbe, 0x09, 0x33, 0xc4,
                    0xa4, 0x1a, 0xc3, 0x51, 0x33, 0x50, 0xee, 0xf3, 0x68, 0x50, 0xb8, 0x04, 0x70,
                    0x8e, 0x1f, 0x33, 0x89, 0x32, 0xce, 0x01, 0xb6, 0x55, 0xa1, 0x63, 0x34, 0x4a,
                    0x45, 0x00, 0xb2, 0x81, 0xc8, 0x75, 0x0c, 0x46, 0x1f, 0x00, 0x38, 0xe7, 0x6f,
                    0x31, 0xb5, 0xae, 0xf9, 0xd7, 0xc8, 0xf1, 0x61, 0x6d, 0x24, 0x46, 0xfb, 0x1f,
                    0x03, 0x80, 0xac, 0xa5, 0x86, 0xb9, 0x26, 0x8a, 0x11, 0x5e, 0x8d, 0x19, 0x1c,
                    0x5e, 0x47, 0xed, 0x04, 0xc4, 0xb7, 0x7b, 0x72, 0x39, 0x47, 0x40, 0x02, 0x57,
                    0x06, 0x84, 0x5c, 0x9c, 0xb7,
                ],
            ),
        ];

        for (k, expected) in COMPRESSED {
            let p = (Point::GENERATOR * &s(*k)).to_affine().unwrap();
            assert_eq!(&p.to_compressed(), expected, "k = {}", k);
            assert_eq!(PointAffine::from_compressed(expected).unwrap(), p);
        }
        for (k, expected) in UNCOMPRESSED {
            let p = (Point::GENERATOR * &s(*k)).to_affine().unwrap();
            assert_eq!(&p.to_uncompressed(), expected, "k = {}", k);
            assert_eq!(PointAffine::from_uncompressed(expected).unwrap(), p);
        }
    }
}
