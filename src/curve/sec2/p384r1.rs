//! Curve p384r1 as defined over the prime field of order 2^384 - 2^128 - 2^96 + 2^32 - 1

use crate::curve::fiat::p384_64::*;
use crate::curve::fiat::p384_scalar_64::*;
use crate::curve::field::{Field, FieldSqrt, Sign};
use crate::curve::{affine, projective, weierstrass::WeierstrassCurve};
use crate::mp::ct::{Choice, CtEqual, CtOption, CtZero};
use crate::params::sec2::p384r1::*;
use crate::{fiat_define_weierstrass_curve, fiat_define_weierstrass_points};
use crate::{fiat_field_ops_impl, fiat_field_sqrt_define};

const GM_LIMBS_SIZE: usize = 6;
const FE_LIMBS_SIZE: usize = 6;

fiat_field_ops_impl!(
    #[doc = "Element of the prime field Fp where p = 2^384 - 2^128 - 2^96 + 2^32 - 1"]
    FieldElement,
    384,
    P_LIMBS,
    FE_LIMBS_SIZE,
    fiat_p384_nonzero,
    fiat_p384_add,
    fiat_p384_sub,
    fiat_p384_mul,
    fiat_p384_square,
    fiat_p384_opp,
    fiat_p384_to_bytes,
    fiat_p384_from_bytes,
    montgomery {
        fiat_p384_to_montgomery,
        fiat_p384_from_montgomery
    }
);
fiat_field_sqrt_define!(FieldElement);

impl FieldElement {
    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    pub fn inverse(&self) -> Self {
        // p-2 = 1*255,0*1,1*32,0*64,1*30,0*1,1*1
        assert!(!self.is_zero());
        let x2 = self.square() * self;
        let x3 = x2.square() * self;
        let x6 = x3.square_rep(3) * &x3;
        let x12 = x6.square_rep(6) * &x6;
        let x15 = x12.square_rep(3) * &x3;
        let x30 = x15.square_rep(15) * &x15;
        let x60 = x30.square_rep(30) * &x30;
        let x120 = x60.square_rep(60) * &x60;
        let x240 = x120.square_rep(120) * &x120;
        let x255 = x240.square_rep(15) * &x15;

        let mut t1 = x255.square_rep(31) * &x30;
        t1 = t1.square_rep(2) * &x2;
        t1 = t1.square_rep(64 + 30) * &x30;
        t1.square_rep(2) * self
    }

    /// Compute the square root 'x' of the field element such that x*x = self
    pub fn sqrt(&self) -> CtOption<Self> {
        // (p+1)/4 = 1*255,0*1,1*32,0*63,1*1,0*30
        let x2 = self.square() * self;
        let x3 = x2.square() * self;
        let x6 = x3.square_rep(3) * &x3;
        let x12 = x6.square_rep(6) * &x6;
        let x15 = x12.square_rep(3) * &x3;
        let x30 = x15.square_rep(15) * &x15;
        let x60 = x30.square_rep(30) * &x30;
        let x120 = x60.square_rep(60) * &x60;
        let x240 = x120.square_rep(120) * &x120;
        let x255 = x240.square_rep(15) * &x15;

        let mut t1 = x255.square_rep(31) * &x30;
        t1 = t1.square_rep(2) * &x2;
        t1 = t1.square_rep(64) * self;
        let r = t1.square_rep(30);

        let r2 = &r * &r;
        CtOption::from((CtEqual::ct_eq(&r2, self), r))
    }
}

fiat_field_ops_impl!(
    #[doc = "Element of the prime field Fp for scalar where p is the order of the SECP384R1 curve"]
    Scalar,
    384,
    ORDER_LIMBS,
    GM_LIMBS_SIZE,
    fiat_p384_scalar_nonzero,
    fiat_p384_scalar_add,
    fiat_p384_scalar_sub,
    fiat_p384_scalar_mul,
    fiat_p384_scalar_square,
    fiat_p384_scalar_opp,
    fiat_p384_scalar_to_bytes,
    fiat_p384_scalar_from_bytes,
    montgomery {
        fiat_p384_scalar_to_montgomery,
        fiat_p384_scalar_from_montgomery
    }
);

impl Scalar {
    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    /// TODO this will change to being a method of NonZeroScalar
    pub fn inverse(&self) -> Self {
        // 1*194,0001110110001101001101100000011111010000110111001011011101111101011000000110100000110110110010010010001011000010100111011110101110110011101100000110010110101011001100110001010010100101110001
        // 1*194,0*3,1*3,0*1,1*2,0*3,1*2,0*1,1*1,0*2,1*2,0*1,1*2,0*6,1*5,0*1,1*1,0*4,1*2,0*1,1*3,0*2,1*1,0*1,1*2,0*1,1*3,0*1,1*5,0*1,1*1,0*1,1*2,0*6,1*2,0*1,1*1,0*5,1*2,0*1,1*2,0*1,1*2,0*2,1*1,0*2,1*1,0*2,1*1,0*3,1*1,0*1,1*2,0*4,1*1,0*1,1*1,0*2,1*3,0*1,1*4,0*1,1*1,0*1,1*3,0*1,1*2,0*2,1*3,0*1,1*2,0*5,1*2,0*2,1*1,0*1,1*2,0*1,1*1,0*1,1*1,0*1,1*2,0*2,1*2,0*2,1*2,0*3,1*1,0*1,1*1,0*2,1*1,0*1,1*1,0*2,1*1,0*1,1*3,0*3,1*1
        assert!(!self.is_zero());
        let b10 = self.square();
        let b11 = &b10 * self; // x2
        let b101 = &b10 * &b11;
        let b111 = &b101 * &b10; // x3
        let b1001 = &b111 * &b10;
        let b1011 = &b1001 * &b10;
        let b1101 = &b1011 * &b10;
        let b1111 = &b1101 * &b10; // x4
        let x4 = &b1111;

        let x8 = x4.square_rep(4) * x4;
        let x16 = x8.square_rep(8) * &x8;
        let x32 = x16.square_rep(16) * &x16;
        let x64 = x32.square_rep(32) * &x32;
        let x96 = x64.square_rep(32) * &x32;
        let x192 = x96.square_rep(96) * &x96;
        let x194 = x192.square_rep(2) * &b11;

        let mut t1 = x194.square_rep(6) * &b111;
        t1 = t1.square_rep(3) * &b11;

        t1 = t1.square_rep(5) * &b11;
        t1 = t1.square_rep(5) * &b1001;
        t1 = t1.square_rep(4) * &b1011;
        t1 = t1.square_rep(10) * &b1111;

        t1 = t1.square_rep(3) * &b101;
        t1 = t1.square_rep(5) * self;
        t1 = t1.square_rep(4) * &b1011;
        t1 = t1.square_rep(4) * &b1001;

        t1 = t1.square_rep(5) * &b1101;
        t1 = t1.square_rep(4) * &b1101;
        t1 = t1.square_rep(4) * &b1111;
        t1 = t1.square_rep(5) * &b1011;

        t1 = t1.square_rep(10) * &b1101;
        t1 = t1.square_rep(9) * &b1101;
        t1 = t1.square_rep(4) * &b1011;
        t1 = t1.square_rep(6) * &b1001;

        t1 = t1.square_rep(3) * self;
        t1 = t1.square_rep(7) * &b1011;
        t1 = t1.square_rep(7) * &b101;
        t1 = t1.square_rep(5) * &b111;

        t1 = t1.square_rep(5) * &b1111;
        t1 = t1.square_rep(5) * &b1011;
        t1 = t1.square_rep(4) * &b1011;
        t1 = t1.square_rep(5) * &b111;

        t1 = t1.square_rep(3) * &b11;
        t1 = t1.square_rep(7) * &b11;
        t1 = t1.square_rep(6) * &b1011;
        t1 = t1.square_rep(4) * &b101;

        t1 = t1.square_rep(3) * &b11;
        t1 = t1.square_rep(4) * &b11;
        t1 = t1.square_rep(4) * &b11;
        t1 = t1.square_rep(6) * &b101;

        t1 = t1.square_rep(5) * &b101;
        t1 = t1.square_rep(5) * &b101;
        t1 = t1.square_rep(2) * &b11;
        t1 = t1.square_rep(4) * self;
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
