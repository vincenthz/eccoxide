//! Curve p192r1 as defined over the prime field of order 2^192 - 2^64 - 1

use crate::curve::fiat::p192r1_64::*;
use crate::curve::fiat::p192r1_scalar_64::*;
use crate::curve::field::{Field, FieldSqrt, Sign};
use crate::curve::{affine, projective, weierstrass::WeierstrassCurve};
use crate::mp::ct::{Choice, CtEqual, CtOption, CtZero};
use crate::params::sec2::p192r1::*;
use crate::{fiat_define_weierstrass_curve, fiat_define_weierstrass_points};
use crate::{fiat_field_ops_impl, fiat_field_sqrt_define};

const GM_LIMBS_SIZE: usize = 3;
const FE_LIMBS_SIZE: usize = 3;

fiat_field_ops_impl!(
    #[doc = "Element of the prime field Fp where p = 2^192 - 2^64 - 1"]
    FieldElement,
    192,
    P_LIMBS,
    FE_LIMBS_SIZE,
    fiat_p192r1_nonzero,
    fiat_p192r1_add,
    fiat_p192r1_sub,
    fiat_p192r1_mul,
    fiat_p192r1_square,
    fiat_p192r1_opp,
    fiat_p192r1_to_bytes,
    fiat_p192r1_from_bytes,
    montgomery {
        fiat_p192r1_to_montgomery,
        fiat_p192r1_from_montgomery
    }
);
fiat_field_sqrt_define!(FieldElement);

impl FieldElement {
    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    pub fn inverse(&self) -> Self {
        assert!(!self.is_zero());
        // 1*127,0*1,1*62,0*1,1*1
        let x2 = self.square() * self;
        let x4 = x2.square_rep(2) * &x2;
        let x6 = x4.square_rep(2) * &x2;
        let x7 = x6.square() * self;
        let x14 = x7.square_rep(7) * &x7;
        let x28 = x14.square_rep(14) * &x14;
        let x56 = x28.square_rep(28) * &x28;
        let x62 = x56.square_rep(6) * &x6;
        let x112 = x56.square_rep(56) * &x56;
        let x126 = x112.square_rep(14) * &x14;
        let x127 = x126.square() * self;

        let mut t1 = x127.square_rep(63) * &x62; // 1*127 0*1 1*62
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1
    }

    /// Compute the square root 'x' of the field element such that x*x = self
    pub fn sqrt(&self) -> CtOption<Self> {
        // 1*128,0*62
        let x2 = self.square() * self;
        let x4 = x2.square_rep(2) * &x2;
        let x8 = x4.square_rep(4) * &x4;
        let x16 = x8.square_rep(8) * &x8;
        let x32 = x16.square_rep(16) * &x16;
        let x64 = x32.square_rep(32) * &x32;
        let x128 = x64.square_rep(64) * &x64;

        let r = x128.square_rep(62); // 1*128 0*62

        let r2 = &r * &r;
        CtOption::from((CtEqual::ct_eq(&r2, self), r))
    }
}

fiat_field_ops_impl!(
    #[doc = "Element of the prime field Fp for scalar where p is the order of the SECP192R1 curve"]
    Scalar,
    192,
    ORDER_LIMBS,
    GM_LIMBS_SIZE,
    fiat_p192r1_scalar_nonzero,
    fiat_p192r1_scalar_add,
    fiat_p192r1_scalar_sub,
    fiat_p192r1_scalar_mul,
    fiat_p192r1_scalar_square,
    fiat_p192r1_scalar_opp,
    fiat_p192r1_scalar_to_bytes,
    fiat_p192r1_scalar_from_bytes,
    montgomery {
        fiat_p192r1_scalar_to_montgomery,
        fiat_p192r1_scalar_from_montgomery
    }
);

impl Scalar {
    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    /// TODO this will change to being a method of NonZeroScalar
    pub fn inverse(&self) -> Self {
        assert!(!self.is_zero());
        let x2 = self.square() * self;
        let x3 = x2.square() * self;
        let x4 = x3.square() * self;
        let x5 = x4.square() * self;
        let x8 = x4.square_rep(4) * &x4;
        let x16 = x8.square_rep(8) * &x8;
        let x32 = x16.square_rep(16) * &x16;
        let x48 = x32.square_rep(16) * &x16;
        let x96 = x48.square_rep(48) * &x48;
        let x97 = x96.square() * self;

        // 1*97,0*2,1*2,0*2,1*3,0*1,1*4,0*1,1*5,0*5,1*2,0*1,1*2,0*4,1*1,0*1,1*1,0*3,1*2,0*1,1*1,0*1,1*4,0*2,1*1,0*2,1*2,0*1,1*2,0*3,1*2,0*1,1*2,0*1,1*1,0*2,1*2,0*1,1*1,0*2,1*1,0*3,1*1,0*1,1*1,0*5,1*1,0*1,1*4
        let mut t1 = x97.square_rep(4) * &x2; // 1*97 0*2 1*2
        t1 = t1.square_rep(5) * &x3; // 0*2 1*3
        t1 = t1.square_rep(5) * &x4; // 0*1 1*4
        t1 = t1.square_rep(6) * &x5; // 0*1 1*5
        t1 = t1.square_rep(7) * &x2; // 0*5 1*2
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2
        t1 = t1.square_rep(5) * self; // 0*4 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(5) * &x2; // 0*3 1*2
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(5) * &x4; // 0*1 1*4
        t1 = t1.square_rep(3) * self; // 0*2 1*1
        t1 = t1.square_rep(4) * &x2; // 0*2 1*2
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2
        t1 = t1.square_rep(5) * &x2; // 0*3 1*2
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(4) * &x2; // 0*2 1*2
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(3) * self; // 0*2 1*1
        t1 = t1.square_rep(4) * self; // 0*3 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(6) * self; // 0*5 1*1
        t1 = t1.square_rep(5) * &x4; // 0*1 1*4

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
    #[cfg(feature = "fast-u64-scalar-mul")]
    fn scale_u64(&self, other: u64) -> Self {
        Point(self.0.scale_u64(other, Curve))
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
