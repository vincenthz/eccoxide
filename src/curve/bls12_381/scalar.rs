//! BLS12-381 scalar field `Scalar` (a.k.a. `Fr`)
//!
//! `Scalar` is the prime field of order `r`, the 255-bit order of the
//! prime-order subgroups G1 (and G2). It is the field of scalars acting on the
//! curve group, so scalar multiplication of a point takes a `Scalar`.
//!
//! As for [`super::fp`], the arithmetic comes from the fiat-crypto Montgomery
//! backend and this module only wires it into the generic [`Field`]
//! abstraction.

use crate::curve::fiat::bls12_381_scalar_64::*;
use crate::curve::field::{Field, Sign};
use crate::fiat_field_montgomery_impl;
use crate::mp::ct::{Choice, CtEqual, CtZero};
use crate::params::bls12_381::ORDER_LIMBS;

/// Number of 64-bit limbs of a scalar field element (255 bits -> 4 limbs).
const GM_LIMBS_SIZE: usize = 4;

fiat_field_montgomery_impl!(
    #[doc = "Element of the BLS12-381 scalar prime field Fr where r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001"]
    Scalar,
    255,
    ORDER_LIMBS,
    GM_LIMBS_SIZE,
    fiat_bls12_381_scalar_non_montgomery_domain_field_element,
    fiat_bls12_381_scalar_nonzero,
    fiat_bls12_381_scalar_add,
    fiat_bls12_381_scalar_sub,
    fiat_bls12_381_scalar_mul,
    fiat_bls12_381_scalar_square,
    fiat_bls12_381_scalar_opp,
    fiat_bls12_381_scalar_to_bytes,
    fiat_bls12_381_scalar_from_bytes,
    fiat_bls12_381_scalar_montgomery_domain_field_element,
    fiat_bls12_381_scalar_to_montgomery,
    fiat_bls12_381_scalar_from_montgomery,
    fiat_bls12_381_scalar_selectznz,
    fiat_bls12_381_scalar_msat,
    fiat_bls12_381_scalar_divstep,
    fiat_bls12_381_scalar_divstep_precomp
);

impl Scalar {
    /// Get the multiplicative inverse.
    ///
    /// This currently delegates to the representation-agnostic Bernstein-Yang
    /// "safegcd" inversion ([`Self::inverse_safegcd`]); see the note on
    /// [`super::fp::Fp::inverse`].
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a
    /// panic.
    pub fn inverse(&self) -> Self {
        self.inverse_safegcd()
    }
}

#[cfg(test)]
mod tests {
    use super::Scalar;
    use crate::fiat_field_unittest;

    fiat_field_unittest!(Scalar);
    crate::fiat_field_safegcd_unittest!(Scalar);
}
