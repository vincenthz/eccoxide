//! Curve p224k1 as defined over the prime field of order  2^224 - 2^32 - 2^12 - 2^11 - 2^9 - 2^7 - 2^4 - 2 - 1

use crate::curve::fiat::p224k1_64::*;
use crate::curve::fiat::p224k1_scalar_64::*;
use crate::curve::field::{Field, FieldSqrt, Sign};
use crate::curve::{
    affine, projective,
    weierstrass::{WeierstrassCurve, WeierstrassCurveA0},
};
use crate::mp::ct::{Choice, CtEqual, CtOption, CtZero};
use crate::params::sec2::p224k1::*;
use crate::{fiat_define_weierstrass_curve, fiat_define_weierstrass_points};
use crate::{fiat_field_ops_impl, fiat_field_sqrt_define};

const GM_LIMBS_SIZE: usize = 4;
const FE_LIMBS_SIZE: usize = 4;

fiat_field_ops_impl!(
    #[doc = "Element of the prime field Fp where p = 2^224 - 2^32 - 2^12 - 2^11 - 2^9 - 2^7 - 2^4 - 2 - 1"]
    FieldElement,
    224,
    P_LIMBS,
    FE_LIMBS_SIZE,
    fiat_p224k1_nonzero,
    fiat_p224k1_add,
    fiat_p224k1_sub,
    fiat_p224k1_mul,
    fiat_p224k1_square,
    fiat_p224k1_opp,
    fiat_p224k1_to_bytes,
    fiat_p224k1_from_bytes,
    montgomery {
        fiat_p224k1_to_montgomery,
        fiat_p224k1_from_montgomery
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
        let x4 = x3.square() * self;
        let x6 = x3.square_rep(3) * &x3;
        let x9 = x6.square_rep(3) * &x3;
        let x11 = x9.square_rep(2) * &x2;
        let x17 = x11.square_rep(6) * &x6;
        let x19 = x17.square_rep(2) * &x2;
        let x22 = x11.square_rep(11) * &x11;
        let x44 = x22.square_rep(22) * &x22;
        let x88 = x44.square_rep(44) * &x44;
        let x176 = x88.square_rep(88) * &x88;
        let x187 = x176.square_rep(11) * &x11;
        let x191 = x187.square_rep(4) * &x4;

        // 1*191,0*1,1*19,0*2,1*1,0*1,1*1,0*1,1*2,0*1,1*1,0*1,1*2
        let mut t1 = x191.square_rep(20) * &x19; // 1*191 0*1 1*19
        t1 = t1.square_rep(3) * self; // 0*2 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2
        t1
    }

    /// Compute the square root 'x' of the field element such that x*x = self
    pub fn sqrt(&self) -> CtOption<Self> {
        // this is using the P == 5 mod 8 method

        // raise 2*self == x to power of (p-5)/8
        let x = self.double();
        let gamma = {
            let x2 = x.square() * &x;
            let x3 = x2.square() * &x;
            let x4 = x3.square() * &x;
            let x6 = x3.square_rep(3) * &x3;
            let x9 = x6.square_rep(3) * &x3;
            let x11 = x9.square_rep(2) * &x2;
            let x17 = x11.square_rep(6) * &x6;
            let x19 = x17.square_rep(2) * &x2;
            let x22 = x11.square_rep(11) * &x11;
            let x44 = x22.square_rep(22) * &x22;
            let x88 = x44.square_rep(44) * &x44;
            let x176 = x88.square_rep(88) * &x88;
            let x187 = x176.square_rep(11) * &x11;
            let x191 = x187.square_rep(4) * &x4;

            // 1*191,0*1,1*19,0*2,1*1,0*1,1*1,0*1,1*2,0*1,1*1
            let mut t1 = x191.square_rep(20) * &x19; // 1*191 0*1 1*19
            t1 = t1.square_rep(3) * &x; // 0*2 1*1
            t1 = t1.square_rep(2) * &x; // 0*1 1*1
            t1 = t1.square_rep(3) * &x2; // 0*1 1*2
            t1 = t1.square_rep(2) * &x; // 0*1 1*1
            t1
        };

        let ggamma = self * &gamma;
        let i = (&ggamma * &gamma).double();
        let im1 = i - Self::one();
        let r = ggamma * im1;

        let r2 = &r * &r;
        CtOption::from((CtEqual::ct_eq(&r2, self), r))
    }
}

fiat_field_ops_impl!(
    #[doc = "Element of the prime field Fp for scalar where p is the order of the SECP224K1 curve"]
    Scalar,
    225,
    ORDER_LIMBS,
    GM_LIMBS_SIZE,
    fiat_p224k1_scalar_nonzero,
    fiat_p224k1_scalar_add,
    fiat_p224k1_scalar_sub,
    fiat_p224k1_scalar_mul,
    fiat_p224k1_scalar_square,
    fiat_p224k1_scalar_opp,
    fiat_p224k1_scalar_to_bytes,
    fiat_p224k1_scalar_from_bytes,
    montgomery {
        fiat_p224k1_scalar_to_montgomery,
        fiat_p224k1_scalar_from_montgomery
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
        let x6 = x5.square() * self;

        // 1*1,0*111,1*3,0*1,1*3,0*2,1*3,0*1,1*1,0*3,1*2,0*1,1*1,0*2,1*1,0*1,1*3,0*1,1*2,0*3,1*2,0*4,1*2,0*4,1*1,0*2,1*2,0*2,1*1,0*1,1*1,0*1,1*4,0*4,1*1,0*1,1*1,0*1,1*1,0*2,1*1,0*1,1*3,0*3,1*1,0*1,1*3,0*1,1*2,0*1,1*1,0*2,1*6,0*1,1*2,0*3,1*5,0*1,1*1,0*1,1*1
        let mut t1 = self.square_rep(114) * &x3; // 1*1 0*111 1*3
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(5) * &x3; // 0*2 1*3
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(5) * &x2; // 0*3 1*2
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(3) * self; // 0*2 1*1
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2
        t1 = t1.square_rep(5) * &x2; // 0*3 1*2
        t1 = t1.square_rep(6) * &x2; // 0*4 1*2
        t1 = t1.square_rep(5) * self; // 0*4 1*1
        t1 = t1.square_rep(4) * &x2; // 0*2 1*2
        t1 = t1.square_rep(3) * self; // 0*2 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(5) * &x4; // 0*1 1*4
        t1 = t1.square_rep(5) * self; // 0*4 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(3) * self; // 0*2 1*1
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(4) * self; // 0*3 1*1
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(8) * &x6; // 0*2 1*6
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2
        t1 = t1.square_rep(8) * &x5; // 0*3 1*5
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
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
    #[cfg(feature = "fast-u64-scalar-mul")]
    fn scale_u64(&self, other: u64) -> Self {
        Point(self.0.scale_a0_u64(other, Curve))
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
