use crate::curve::fiat::secp256k1_64::*;
use crate::curve::fiat::secp256k1_scalar_64::*;
use crate::curve::field::{Field, FieldSqrt, Sign};
use crate::curve::{
    affine, projective,
    weierstrass::{WeierstrassCurve, WeierstrassCurveA0},
};
use crate::mp::ct::{Choice, CtEqual, CtOption, CtZero};
use crate::params::sec2::p256k1::*;
use crate::{fiat_define_weierstrass_curve, fiat_define_weierstrass_points};
use crate::{fiat_field_ops_impl, fiat_field_sqrt_define};

const GM_LIMBS_SIZE: usize = 4;
const FE_LIMBS_SIZE: usize = 4;

fiat_field_ops_impl!(
    FieldElement,
    256,
    P_LIMBS,
    FE_LIMBS_SIZE,
    fiat_secp256k1_nonzero,
    fiat_secp256k1_add,
    fiat_secp256k1_sub,
    fiat_secp256k1_mul,
    fiat_secp256k1_square,
    fiat_secp256k1_opp,
    fiat_secp256k1_to_bytes,
    fiat_secp256k1_from_bytes,
    montgomery {
        fiat_secp256k1_to_montgomery,
        fiat_secp256k1_from_montgomery
    }
);
fiat_field_sqrt_define!(FieldElement);

impl FieldElement {
    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    pub fn inverse(&self) -> Self {
        assert!(!self.is_zero());
        let x2 = self.square() * self;
        let x3 = x2.square() * self;
        let x6 = x3.square_rep(3) * &x3;
        let x9 = x6.square_rep(3) * &x3;
        let x11 = x9.square_rep(2) * &x2;
        let x22 = x11.square_rep(11) * &x11;
        let x44 = x22.square_rep(22) * &x22;
        let x88 = x44.square_rep(44) * &x44;
        let x176 = x88.square_rep(88) * &x88;
        let x220 = x176.square_rep(44) * &x44;
        let x223 = x220.square_rep(3) * &x3;

        let mut t1 = x223.square_rep(23) * &x22;
        t1 = t1.square_rep(5) * self;
        t1 = t1.square_rep(3) * &x2;
        t1 = t1.square_rep(2);

        t1 * self
    }

    /// Compute the square root 'x' of the field element such that x*x = self
    pub fn sqrt(&self) -> CtOption<Self> {
        let x2 = self.square() * self;
        let x3 = x2.square() * self;
        let x6 = x3.square_rep(3) * &x3;
        let x9 = x6.square_rep(3) * &x3;
        let x11 = x9.square_rep(2) * &x2;
        let x22 = x11.square_rep(11) * &x11;
        let x44 = x22.square_rep(22) * &x22;
        let x88 = x44.square_rep(44) * &x44;
        let x176 = x88.square_rep(88) * &x88;
        let x220 = x176.square_rep(44) * &x44;
        let x223 = x220.square_rep(3) * &x3;

        let mut t1 = x223.square_rep(23) * &x22;
        t1 = t1.square_rep(6) * &x2;
        t1 = &t1 * &t1;

        let r = &t1 * &t1;
        let r2 = &r * &r;
        CtOption::from((CtEqual::ct_eq(&r2, self), r))
    }
}

fiat_field_ops_impl!(
    Scalar,
    256,
    ORDER_LIMBS,
    GM_LIMBS_SIZE,
    fiat_secp256k1_scalar_nonzero,
    fiat_secp256k1_scalar_add,
    fiat_secp256k1_scalar_sub,
    fiat_secp256k1_scalar_mul,
    fiat_secp256k1_scalar_square,
    fiat_secp256k1_scalar_opp,
    fiat_secp256k1_scalar_to_bytes,
    fiat_secp256k1_scalar_from_bytes,
    montgomery {
        fiat_secp256k1_scalar_to_montgomery,
        fiat_secp256k1_scalar_from_montgomery
    }
);

impl Scalar {
    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    /// TODO this will change to being a method of NonZeroScalar
    pub fn inverse(&self) -> Self {
        assert!(!self.is_zero());
        let x = self;
        let u2 = self.square();
        let x2 = &u2 * x;
        let u5 = &u2 * &x2;
        let x3 = &u5 * &u2;
        let u9 = &x3 * &u2;
        let u11 = &u9 * &u2;
        let u13 = &u11 * &u2;

        let x6 = u13.square().square() * &u11;
        let x8 = x6.square().square() * &x2;

        let x14 = x8.square_rep(6) * &x6;
        let x28 = x14.square_rep(14) * &x14;
        let x56 = x28.square_rep(28) * &x28;
        let x112 = x56.square_rep(56) * &x56;
        let x126 = x112.square_rep(14) * &x14;

        let mut t = x126.square_rep(3) * &u5;
        t = t.square_rep(4) * &x3;
        t = t.square_rep(4) * &u5;
        t = t.square_rep(5) * &u11;
        t = t.square_rep(4) * &u11;
        t = t.square_rep(4) * &x3;
        t = t.square_rep(5) * &x3;
        t = t.square_rep(6) * &u13;
        t = t.square_rep(4) * &u5;
        t = t.square_rep(3) * &x3;
        t = t.square_rep(5) * &u9;

        t = t.square_rep(6) * &u5;
        t = t.square_rep(10) * &x3;
        t = t.square_rep(4) * &x3;
        t = t.square_rep(9) * &x8;
        t = t.square_rep(5) * &u9;
        t = t.square_rep(6) * &u11;
        t = t.square_rep(4) * &u13;
        t = t.square_rep(5) * &x2;
        t = t.square_rep(6) * &u13;
        t = t.square_rep(10) * &u13;
        t = t.square_rep(4) * &u9;
        t = t.square_rep(6) * x;
        t = t.square_rep(8) * &x6;

        t
    }
}

fiat_define_weierstrass_curve!(FieldElement);
fiat_define_weierstrass_points!(FieldElement);

impl WeierstrassCurveA0 for Curve {}

impl Point {
    fn add_or_double<'b>(&self, other: &'b Point) -> Point {
        Point(self.0.add_or_double_a0(&other.0, Curve))
    }
    fn scale<'b>(&self, other: &'b Scalar) -> Self {
        Point(self.0.scale_a0(&other.to_bytes(), Curve))
    }
}

#[cfg(test)]
mod tests {
    mod fe {
        use super::super::FieldElement;
        use crate::{fiat_field_sqrt_unittest, fiat_field_unittest};

        fiat_field_unittest!(FieldElement);
        fiat_field_sqrt_unittest!(FieldElement);
    }
    mod gm {
        use super::super::Scalar;
        use crate::fiat_field_unittest;
        fiat_field_unittest!(Scalar);
    }
}
