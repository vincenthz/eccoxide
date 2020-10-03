//! Curve p192k1 as defined over the prime field of order 2^192 - 2^32 - 2^12 - 2^8 - 2^7 - 2^6 - 2^3 - 1

use crate::curve::fiat::p192k1_64::*;
use crate::curve::fiat::p192k1_scalar_64::*;
use crate::curve::field::{Field, FieldSqrt, Sign};
use crate::curve::{
    affine, projective,
    weierstrass::{WeierstrassCurve, WeierstrassCurveA0},
};
use crate::mp::ct::{Choice, CtEqual, CtOption, CtZero};
use crate::params::sec2::p192k1::*;
use crate::{fiat_define_weierstrass_curve, fiat_define_weierstrass_points};
use crate::{fiat_field_ops_impl, fiat_field_sqrt_define};

const GM_LIMBS_SIZE: usize = 3;
const FE_LIMBS_SIZE: usize = 3;

fiat_field_ops_impl!(
    #[doc = "Element of the prime field Fp where p = 2^192 - 2^32 - 2^12 - 2^8 - 2^7 - 2^6 - 2^3 - 1"]
    FieldElement,
    192,
    P_LIMBS,
    FE_LIMBS_SIZE,
    fiat_p192k1_nonzero,
    fiat_p192k1_add,
    fiat_p192k1_sub,
    fiat_p192k1_mul,
    fiat_p192k1_square,
    fiat_p192k1_opp,
    fiat_p192k1_to_bytes,
    fiat_p192k1_from_bytes,
    montgomery {
        fiat_p192k1_to_montgomery,
        fiat_p192k1_from_montgomery
    }
);
fiat_field_sqrt_define!(FieldElement);

impl FieldElement {
    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    pub fn inverse(&self) -> Self {
        // 1*159,0*1,1*19,0*1,1*3,0*3,1*2,0*1,1*1,0*1,1*1
        assert!(!self.is_zero());
        let x2 = self.square() * self;
        let x3 = x2.square() * self;
        let x6 = x3.square_rep(3) * &x3;
        let x9 = x6.square_rep(3) * &x3;
        let x18 = x9.square_rep(9) * &x9;
        let x19 = x18.square() * self;
        let x38 = x19.square_rep(19) * &x19;
        let x76 = x38.square_rep(38) * &x38;
        let x152 = x76.square_rep(76) * &x76;
        let x158 = x152.square_rep(6) * &x6;
        let x159 = x158.square() * self;
        let mut t1 = x159.square_rep(20) * &x19; // 1*159 0*1 1*19
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(5) * &x2; // 0*3 1*2
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1
    }

    /// Compute the square root 'x' of the field element such that x*x = self
    pub fn sqrt(&self) -> CtOption<Self> {
        // 1*159,0*1,1*19,0*1,1*3,0*3,1*3,0*1
        let x2 = self.square() * self;
        let x3 = x2.square() * self;
        let x6 = x3.square_rep(3) * &x3;
        let x9 = x6.square_rep(3) * &x3;
        let x18 = x9.square_rep(9) * &x9;
        let x19 = x18.square() * self;
        let x38 = x19.square_rep(19) * &x19;
        let x76 = x38.square_rep(38) * &x38;
        let x152 = x76.square_rep(76) * &x76;
        let x158 = x152.square_rep(6) * &x6;
        let x159 = x158.square() * self;

        let mut t1 = x159.square_rep(20) * &x19; // 1*159 0*1 1*19
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(6) * &x3; // 0*3 1*3

        let r = &t1 * &t1;
        let r2 = &r * &r;
        CtOption::from((CtEqual::ct_eq(&r2, self), r))
    }
}

fiat_field_ops_impl!(
    #[doc = "Element of the prime field Fp for scalar where p is the order of the SECP192K1 curve"]
    Scalar,
    192,
    ORDER_LIMBS,
    GM_LIMBS_SIZE,
    fiat_p192k1_scalar_nonzero,
    fiat_p192k1_scalar_add,
    fiat_p192k1_scalar_sub,
    fiat_p192k1_scalar_mul,
    fiat_p192k1_scalar_square,
    fiat_p192k1_scalar_opp,
    fiat_p192k1_scalar_to_bytes,
    fiat_p192k1_scalar_from_bytes,
    montgomery {
        fiat_p192k1_scalar_to_montgomery,
        fiat_p192k1_scalar_from_montgomery
    }
);

impl Scalar {
    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    /// TODO this will change to being a method of NonZeroScalar
    pub fn inverse(&self) -> Self {
        // 1*95,0*3,1*1,0*2,1*2,0*1,1*4,0*2,1*1,0*1,1*6,0*5,1*1,0*1,1*3,0*4,1*4,0*1,1*2,0*1,1*1,0*2,1*1,0*1,1*1,0*3,1*2,0*2,1*2,0*1,1*1,0*1,1*1,0*2,1*3,0*1,1*1,0*2,1*2,0*1,1*4,0*1,1*6,0*1,1*2,0*3,1*1,0*1,1*2
        assert!(!self.is_zero());
        let x2 = self.square() * self;
        let x3 = x2.square() * self;
        let x4 = x3.square() * self;
        let x6 = x3.square_rep(3) * &x3;
        let x9 = x6.square_rep(3) * &x3;
        let x18 = x9.square_rep(9) * &x9;
        let x19 = x18.square() * self;
        let x38 = x19.square_rep(19) * &x19;
        let x76 = x38.square_rep(38) * &x38;
        let x95 = x76.square_rep(19) * &x19;

        let mut t1 = x95.square_rep(4) * self; // 1*95 0*3 1*1
        t1 = t1.square_rep(4) * &x2; // 0*2 1*2
        t1 = t1.square_rep(5) * &x4; // 0*1 1*4
        t1 = t1.square_rep(3) * self; // 0*2 1*1
        t1 = t1.square_rep(7) * &x6; // 0*1 1*6
        t1 = t1.square_rep(6) * self; // 0*5 1*1
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(8) * &x4; // 0*4 1*4
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(3) * self; // 0*2 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(5) * &x2; // 0*3 1*2
        t1 = t1.square_rep(4) * &x2; // 0*2 1*2
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(5) * &x3; // 0*2 1*3
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(4) * &x2; // 0*2 1*2
        t1 = t1.square_rep(5) * &x4; // 0*1 1*4
        t1 = t1.square_rep(7) * &x6; // 0*1 1*6
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2
        t1 = t1.square_rep(4) * self; // 0*3 1*1
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2

        t1
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
