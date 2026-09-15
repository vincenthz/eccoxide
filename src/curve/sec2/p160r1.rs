//! Curve p160r1 as defined in SEC2.

use crate::curve::fiat::p160r1_64::*;
use crate::curve::fiat::p160r1_scalar_64::*;
use crate::curve::field::{Field, FieldSqrt, Sign};
use crate::curve::{
    affine, projective,
    weierstrass::{WeierstrassCurve, WeierstrassCurveAM3},
};
use crate::mp::ct::{Choice, CtEqual, CtOption, CtZero};
use crate::params::sec2::p160r1::*;
use crate::{
    fiat_define_weierstrass_curve, fiat_define_weierstrass_curve_am3,
    fiat_define_weierstrass_points,
};
use crate::{fiat_field_montgomery_impl, fiat_field_sqrt_define};

const GM_LIMBS_SIZE: usize = 3;
const FE_LIMBS_SIZE: usize = 3;

fiat_field_montgomery_impl!(
    #[doc = "Element of the 160-bit prime field of the SECP160R1 curve"]
    FieldElement,
    160,
    P_LIMBS,
    FE_LIMBS_SIZE,
    fiat_p160r1_non_montgomery_domain_field_element,
    fiat_p160r1_nonzero,
    fiat_p160r1_add,
    fiat_p160r1_sub,
    fiat_p160r1_mul,
    fiat_p160r1_square,
    fiat_p160r1_opp,
    fiat_p160r1_to_bytes,
    fiat_p160r1_from_bytes,
    fiat_p160r1_montgomery_domain_field_element,
    fiat_p160r1_to_montgomery,
    fiat_p160r1_from_montgomery,
    fiat_p160r1_selectznz,
    fiat_p160r1_msat,
    fiat_p160r1_divstep,
    fiat_p160r1_divstep_precomp
);
fiat_field_sqrt_define!(FieldElement);

impl FieldElement {
    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    pub fn inverse(&self) -> Self {
        self.inverse_safegcd()
    }

    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    pub fn inverse_fermat(&self) -> Self {
        assert!(!self.is_zero());
        // Fixed addition chain for modulus - 2.
        let x2 = self.square() * self;
        let x4 = x2.square_rep(2) * &x2;
        let x8 = x4.square_rep(4) * &x4;
        let x16 = x8.square_rep(8) * &x8;
        let x32 = x16.square_rep(16) * &x16;
        let x64 = x32.square_rep(32) * &x32;
        let x128 = x64.square_rep(64) * &x64;
        let x3 = x2.square() * self;
        let x7 = x4.square_rep(3) * &x3;
        let x15 = x8.square_rep(7) * &x7;
        let x14 = x7.square_rep(7) * &x7;
        let x29 = x15.square_rep(14) * &x14;
        let mut t1 = x128.clone();
        t1 = t1.square_rep(30) * &x29;
        t1 = t1.square_rep(2) * self;
        t1
    }

    /// Compute the square root 'x' of the field element such that x*x = self
    pub fn sqrt(&self) -> CtOption<Self> {
        // p = 3 mod 4: raise to (p + 1) / 4, then check the result.
        let x2 = self.square() * self;
        let x3 = x2.square() * self;
        let x5 = x3.square_rep(2) * &x2;
        let x4 = x2.square_rep(2) * &x2;
        let x9 = x5.square_rep(4) * &x4;
        let x8 = x4.square_rep(4) * &x4;
        let x17 = x9.square_rep(8) * &x8;
        let x16 = x8.square_rep(8) * &x8;
        let x33 = x17.square_rep(16) * &x16;
        let x32 = x16.square_rep(16) * &x16;
        let x65 = x33.square_rep(32) * &x32;
        let x64 = x32.square_rep(32) * &x32;
        let x129 = x65.square_rep(64) * &x64;
        let mut t1 = x129.clone();
        t1 = t1.square_rep(29);

        let r2 = t1.square();
        CtOption::from((CtEqual::ct_eq(&r2, self), t1))
    }
}

fiat_field_montgomery_impl!(
    #[doc = "Element of the prime field Fp for scalar where p is the order of the SECP160R1 curve"]
    Scalar,
    161,
    ORDER_LIMBS,
    GM_LIMBS_SIZE,
    fiat_p160r1_scalar_non_montgomery_domain_field_element,
    fiat_p160r1_scalar_nonzero,
    fiat_p160r1_scalar_add,
    fiat_p160r1_scalar_sub,
    fiat_p160r1_scalar_mul,
    fiat_p160r1_scalar_square,
    fiat_p160r1_scalar_opp,
    fiat_p160r1_scalar_to_bytes,
    fiat_p160r1_scalar_from_bytes,
    fiat_p160r1_scalar_montgomery_domain_field_element,
    fiat_p160r1_scalar_to_montgomery,
    fiat_p160r1_scalar_from_montgomery,
    fiat_p160r1_scalar_selectznz,
    fiat_p160r1_scalar_msat,
    fiat_p160r1_scalar_divstep,
    fiat_p160r1_scalar_divstep_precomp
);

impl Scalar {
    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    pub fn inverse(&self) -> Self {
        self.inverse_safegcd()
    }

    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    /// TODO this will change to being a method of NonZeroScalar
    pub fn inverse_fermat(&self) -> Self {
        assert!(!self.is_zero());
        // Fixed addition chain for modulus - 2.
        let x2 = self.square() * self;
        let x3 = x2.square() * self;
        let x5 = x3.square_rep(2) * &x2;
        let x4 = x2.square_rep(2) * &x2;
        let mut t1 = self.clone();
        t1 = t1.square_rep(84) * &x5;
        t1 = t1.square_rep(2) * self;
        t1 = t1.square_rep(4) * &x2;
        t1 = t1.square_rep(3) * self;
        t1 = t1.square_rep(8) * &x5;
        t1 = t1.square_rep(3) * self;
        t1 = t1.square_rep(3) * self;
        t1 = t1.square_rep(6) * &x4;
        t1 = t1.square_rep(2) * self;
        t1 = t1.square_rep(4) * &x3;
        t1 = t1.square_rep(3) * &x2;
        t1 = t1.square_rep(2) * self;
        t1 = t1.square_rep(6) * &x4;
        t1 = t1.square_rep(3) * self;
        t1 = t1.square_rep(2) * self;
        t1 = t1.square_rep(5) * &x3;
        t1 = t1.square_rep(2) * self;
        t1 = t1.square_rep(2) * self;
        t1 = t1.square_rep(3) * self;
        t1 = t1.square_rep(4) * self;
        t1 = t1.square_rep(3) * self;
        t1 = t1.square_rep(2) * self;
        t1 = t1.square_rep(2) * self;
        t1 = t1.square_rep(2) * self;
        t1
    }
}

fiat_define_weierstrass_curve!(FieldElement);
fiat_define_weierstrass_points!(FieldElement);
fiat_define_weierstrass_curve_am3!(FieldElement);

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

    mod point {
        crate::fiat_curve_point_unittest!();
    }
}
