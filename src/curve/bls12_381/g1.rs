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
    fiat_define_weierstrass_curve, fiat_define_weierstrass_curve_a0, fiat_define_weierstrass_points,
};

fiat_define_weierstrass_curve!(FieldElement);
fiat_define_weierstrass_curve_a0!(FieldElement);
fiat_define_weierstrass_points!(FieldElement);

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
        let recovered = PointAffine::decompress(&x, sign).unwrap();
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
}
