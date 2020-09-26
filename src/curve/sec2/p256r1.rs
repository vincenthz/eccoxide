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
    FieldElement,
    256,
    P_LIMBS,
    FE_LIMBS_SIZE,
    fiat_p256_nonzero,
    fiat_p256_to_montgomery,
    fiat_p256_from_montgomery,
    fiat_p256_add,
    fiat_p256_sub,
    fiat_p256_mul,
    fiat_p256_square,
    fiat_p256_opp,
    fiat_p256_to_bytes,
    fiat_p256_from_bytes
);
fiat_field_sqrt_define!(FieldElement);

impl FieldElement {
    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    pub fn inverse(&self) -> Self {
        assert!(!self.is_zero());
        todo!()
    }

    /// Compute the square root 'x' of the field element such that x*x = self
    pub fn sqrt(&self) -> CtOption<Self> {
        todo!()
    }
}

fiat_field_ops_impl!(
    Scalar,
    256,
    ORDER_LIMBS,
    GM_LIMBS_SIZE,
    fiat_p256_scalar_nonzero,
    fiat_p256_scalar_to_montgomery,
    fiat_p256_scalar_from_montgomery,
    fiat_p256_scalar_add,
    fiat_p256_scalar_sub,
    fiat_p256_scalar_mul,
    fiat_p256_scalar_square,
    fiat_p256_scalar_opp,
    fiat_p256_scalar_to_bytes,
    fiat_p256_scalar_from_bytes
);

impl Scalar {
    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    /// TODO this will change to being a method of NonZeroScalar
    pub fn inverse(&self) -> Self {
        assert!(!self.is_zero());
        todo!()
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
