//! Curve p521r1 as defined over the prime field of order 2^521 - 1

use crate::curve::fiat::p521_64::*;
use crate::curve::fiat::p521_scalar_64::*;
use crate::curve::field::{Field, FieldSqrt, Sign};
use crate::curve::{affine, projective, weierstrass::WeierstrassCurve};
use crate::mp::ct::{Choice, CtEqual, CtOption, CtZero};
use crate::params::sec2::p521r1::*;
use crate::{fiat_define_weierstrass_curve, fiat_define_weierstrass_points};
use crate::{fiat_field_ops_impl, fiat_field_sqrt_define};

const GM_LIMBS_SIZE: usize = 9;
const FE_LIMBS_SIZE: usize = 9;

fiat_field_ops_impl!(
    #[doc = "Element of the prime field Fp where p = 2^521-1"]
    FieldElement,
    521,
    P_LIMBS,
    FE_LIMBS_SIZE,
    fiat_p521_nonzero,
    fiat_p521_add,
    fiat_p521_sub,
    fiat_p521_mul,
    fiat_p521_square,
    fiat_p521_opp,
    fiat_p521_to_bytes,
    fiat_p521_from_bytes,
    montgomery {
        fiat_p521_to_montgomery,
        fiat_p521_from_montgomery
    }
);
fiat_field_sqrt_define!(FieldElement);

impl FieldElement {
    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    pub fn inverse(&self) -> Self {
        // p-2 = 1*519,0*1,1*1
        assert!(!self.is_zero());

        let x2 = self.square() * self;
        let x3 = x2.square() * self;
        let x6 = x3.square_rep(3) * &x3;
        let x9 = x6.square_rep(3) * &x3;
        let x12 = x6.square_rep(6) * &x6;
        let x24 = x12.square_rep(12) * &x12;
        let x48 = x24.square_rep(24) * &x24;
        let x51 = x48.square_rep(3) * &x3;
        let x102 = x51.square_rep(51) * &x51;
        let x204 = x102.square_rep(102) * &x102;
        let x408 = x204.square_rep(204) * &x204;
        let x510 = x408.square_rep(102) * &x102;
        let x519 = x510.square_rep(9) * &x9;
        x519.square_rep(2) * self
    }

    /// Compute the square root 'x' of the field element such that x*x = self
    pub fn sqrt(&self) -> CtOption<Self> {
        // (p+1)/4 = 1*1,0*519
        let r = self.square_rep(519);
        let r2 = &r * &r;
        CtOption::from((CtEqual::ct_eq(&r2, self), r))
    }
}

fiat_field_ops_impl!(
    #[doc = "Element of the prime field Fp for scalar where p is the order of the SECP521R1 curve"]
    Scalar,
    521,
    ORDER_LIMBS,
    GM_LIMBS_SIZE,
    fiat_p521_scalar_nonzero,
    fiat_p521_scalar_add,
    fiat_p521_scalar_sub,
    fiat_p521_scalar_mul,
    fiat_p521_scalar_square,
    fiat_p521_scalar_opp,
    fiat_p521_scalar_to_bytes,
    fiat_p521_scalar_from_bytes,
    montgomery {
        fiat_p521_scalar_to_montgomery,
        fiat_p521_scalar_from_montgomery
    }
);

impl Scalar {
    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    /// TODO this will change to being a method of NonZeroScalar
    pub fn inverse(&self) -> Self {
        // p-2 = 1*262,0*1,1*1,0*2,1*1,0*1,1*1,0*3,1*2,0*4,1*2,0*1,1*1,0*4,1*4,0*5,1*3,0*1,1*6,0*2,1*1,0*1,1*5,0*2,1*1,0*1,1*2,0*2,1*2,0*1,1*1,0*1,1*2,0*1,1*9,0*2,1*2,0*9,1*1,0*1,1*1,0*2,1*1,0*3,1*4,0*1,1*3,0*4,1*1,0*2,1*2,0*1,1*1,0*2,1*1,0*1,1*3,0*1,1*1,0*6,1*3,0*1,1*3,0*1,1*2,0*1,1*1,0*1,1*3,0*2,1*1,0*2,1*2,0*1,1*3,0*3,1*1,0*3,1*1,0*2,1*2,0*2,1*3,0*3,1*1,0*3,1*4,0*1,1*1,0*1,1*3,0*1,1*1,0*1,1*3,0*1,1*2,0*1,1*2,0*1,1*5,0*1,1*2,0*1,1*3,0*3,1*4,0*1,1*1,0*2,1*1,0*3,1*1,0*2,1*3,0*4,1*2,0*2,1*1,0*7,1*3
        assert!(!self.is_zero());

        let x2 = self.square() * self;
        let x3 = x2.square() * self;
        let x4 = x3.square() * self;
        let x5 = x4.square() * self;
        let x6 = x5.square() * self;
        let x9 = x6.square_rep(3) * &x3;
        let x15 = x9.square_rep(6) * &x6;
        let x30 = x15.square_rep(15) * &x15;
        let x60 = x30.square_rep(30) * &x30;
        let x65 = x60.square_rep(5) * &x5;
        let x130 = x65.square_rep(65) * &x65;
        let x260 = x130.square_rep(130) * &x130;
        let x262 = x260.square_rep(2) * &x2;

        let mut t1 = x262.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(3) * self; // 0*2 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(5) * &x2; // 0*3 1*2
        t1 = t1.square_rep(6) * &x2; // 0*4 1*2
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(8) * &x4; // 0*4 1*4
        t1 = t1.square_rep(8) * &x3; // 0*5 1*3
        t1 = t1.square_rep(7) * &x6; // 0*1 1*6
        t1 = t1.square_rep(3) * self; // 0*2 1*1
        t1 = t1.square_rep(6) * &x5; // 0*1 1*5
        t1 = t1.square_rep(3) * self; // 0*2 1*1
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2
        t1 = t1.square_rep(4) * &x2; // 0*2 1*2
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2
        t1 = t1.square_rep(10) * &x9; // 0*1 1*9
        t1 = t1.square_rep(4) * &x2; // 0*2 1*2
        t1 = t1.square_rep(10) * self; // 0*9 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(3) * self; // 0*2 1*1
        t1 = t1.square_rep(7) * &x4; // 0*3 1*4
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(5) * self; // 0*4 1*1
        t1 = t1.square_rep(4) * &x2; // 0*2 1*2
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(3) * self; // 0*2 1*1
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(9) * &x3; // 0*6 1*3
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(3) * self; // 0*2 1*1
        t1 = t1.square_rep(4) * &x2; // 0*2 1*2
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(4) * self; // 0*3 1*1
        t1 = t1.square_rep(4) * self; // 0*3 1*1
        t1 = t1.square_rep(4) * &x2; // 0*2 1*2
        t1 = t1.square_rep(5) * &x3; // 0*2 1*3
        t1 = t1.square_rep(4) * self; // 0*3 1*1
        t1 = t1.square_rep(7) * &x4; // 0*3 1*4
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2
        t1 = t1.square_rep(6) * &x5; // 0*1 1*5
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(7) * &x4; // 0*3 1*4
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(3) * self; // 0*2 1*1
        t1 = t1.square_rep(4) * self; // 0*3 1*1
        t1 = t1.square_rep(5) * &x3; // 0*2 1*3
        t1 = t1.square_rep(6) * &x2; // 0*4 1*2
        t1 = t1.square_rep(3) * self; // 0*2 1*1
        t1 = t1.square_rep(10) * &x3; // 0*7 1*3
        t1
    }
}

fiat_define_weierstrass_curve!(FieldElement);
fiat_define_weierstrass_points!(FieldElement);

impl Point {
    fn add_or_double<'b>(&self, other: &'b Point) -> Point {
        Point(self.0.add_or_double(&other.0, Curve))
    }
    fn scale<'b>(&self, other: &'b Scalar) -> Self {
        Point(self.0.scale(&other.to_bytes(), Curve))
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
