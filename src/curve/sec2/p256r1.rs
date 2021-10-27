//! Curve p256r1 as defined over the prime field of order 2^256 - 2^224 + 2^192 + 2^96 - 1
use crate::curve::fiat::p256_64::*;
use crate::curve::fiat::p256_scalar_64::*;
use crate::curve::field::{Field, FieldSqrt, Sign};
use crate::curve::{affine, projective, weierstrass::WeierstrassCurve};
use crate::mp::ct::{Choice, CtEqual, CtOption, CtZero};
use crate::params::sec2::p256r1::*;
use crate::{fiat_define_weierstrass_curve, fiat_define_weierstrass_points};
use crate::{fiat_field_ops_impl, fiat_field_sqrt_define};

const GM_LIMBS_SIZE: usize = 4;
const FE_LIMBS_SIZE: usize = 4;

fiat_field_ops_impl!(
    #[doc = "Element of the prime field Fp where p = 2^256 - 2^224 + 2^192 + 2^96 - 1"]
    FieldElement,
    256,
    P_LIMBS,
    FE_LIMBS_SIZE,
    fiat_p256_nonzero,
    fiat_p256_add,
    fiat_p256_sub,
    fiat_p256_mul,
    fiat_p256_square,
    fiat_p256_opp,
    fiat_p256_to_bytes,
    fiat_p256_from_bytes,
    montgomery {
        fiat_p256_to_montgomery,
        fiat_p256_from_montgomery
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
        let x12 = x6.square_rep(6) * &x6;
        let x15 = x12.square_rep(3) * &x3;
        let x30 = x15.square_rep(15) * &x15;
        let x32 = x30.square_rep(2) * &x2;

        let mut t1 = x32.square_rep(32) * self;
        t1 = t1.square_rep(96 + 32) * &x32;
        t1 = t1.square_rep(32) * &x32;
        t1 = t1.square_rep(30) * &x30;
        t1 = t1.square_rep(2);
        t1 * self
    }

    /// Compute the square root 'x' of the field element such that x*x = self
    pub fn sqrt(&self) -> CtOption<Self> {
        // (p+1)/4 = 1*32,0*31,1*1,0*95,1*1,0*94

        let x2 = self.square() * self;
        let x3 = x2.square() * self;
        let x6 = x3.square_rep(3) * &x3;
        let x12 = x6.square_rep(6) * &x6;
        let x15 = x12.square_rep(3) * &x3;
        let x30 = x15.square_rep(15) * &x15;
        let x32 = x30.square_rep(2) * &x2;

        let mut t1 = x32.square_rep(32) * self;
        t1 = t1.square_rep(96) * self;
        let r = t1.square_rep(94);

        let r2 = &r * &r;
        CtOption::from((CtEqual::ct_eq(&r2, self), r))
    }
}

fiat_field_ops_impl!(
    #[doc = "Element of the prime field Fp for scalar where p is the order of the SECP256R1 curve"]
    Scalar,
    256,
    ORDER_LIMBS,
    GM_LIMBS_SIZE,
    fiat_p256_scalar_nonzero,
    fiat_p256_scalar_add,
    fiat_p256_scalar_sub,
    fiat_p256_scalar_mul,
    fiat_p256_scalar_square,
    fiat_p256_scalar_opp,
    fiat_p256_scalar_to_bytes,
    fiat_p256_scalar_from_bytes,
    montgomery {
        fiat_p256_scalar_to_montgomery,
        fiat_p256_scalar_from_montgomery
    }
);

impl Scalar {
    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    /// TODO this will change to being a method of NonZeroScalar
    pub fn inverse(&self) -> Self {
        assert!(!self.is_zero());
        let b10 = self.square();
        let b11 = &b10 * self; // x2
        let b101 = &b10 * &b11;
        let b111 = &b101 * &b10; // x3

        let b1010 = b101.square();
        let b1111 = &b1010 * &b101; // x4
        let b10101 = b1010.square() * self;
        let b101010 = b10101.square();

        let b101111 = &b101010 * &b101;

        let x6 = &b101010 * &b10101;
        let x8 = x6.square_rep(2) * &b11;
        let x16 = x8.square_rep(8) * &x8;
        let x32 = x16.square_rep(16) * &x16;

        let mut t1 = x32.square_rep(64) * &x32;
        t1 = t1.square_rep(32) * &x32;
        t1 = t1.square_rep(6) * &b101111;
        t1 = t1.square_rep(5) * &b111;

        t1 = t1.square_rep(4) * &b11;
        t1 = t1.square_rep(5) * &b1111;
        t1 = t1.square_rep(5) * &b10101;
        t1 = t1.square_rep(4) * &b101;

        t1 = t1.square_rep(3) * &b101;
        t1 = t1.square_rep(3) * &b101;
        t1 = t1.square_rep(5) * &b111;
        t1 = t1.square_rep(9) * &b101111;

        t1 = t1.square_rep(6) * &b1111;
        t1 = t1.square_rep(2) * self;
        t1 = t1.square_rep(5) * self;
        t1 = t1.square_rep(6) * &b1111;

        t1 = t1.square_rep(5) * &b111;
        t1 = t1.square_rep(4) * &b111;
        t1 = t1.square_rep(5) * &b111;
        t1 = t1.square_rep(5) * &b101;

        t1 = t1.square_rep(3) * &b11;
        t1 = t1.square_rep(10) * &b101111;
        t1 = t1.square_rep(2) * &b11;
        t1 = t1.square_rep(5) * &b11;

        t1 = t1.square_rep(5) * &b11;
        t1 = t1.square_rep(3) * self;
        t1 = t1.square_rep(7) * &b10101;
        t1 = t1.square_rep(6) * &b1111;

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
