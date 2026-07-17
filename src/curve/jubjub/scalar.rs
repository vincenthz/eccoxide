//! Jubjub scalar field `Scalar`
//!
//! `Scalar` is the prime field of order `r`, the 252-bit order of the Jubjub
//! prime-order subgroup. It is the field of scalars acting on the curve group,
//! so scalar multiplication of a [`Point`](super::Point) takes a `Scalar`.
//!
//! As for the other fiat-backed fields, the arithmetic comes from the
//! fiat-crypto Montgomery backend and this module only wires it into the generic
//! [`Field`] abstraction.
//!
//! # Generating the backend
//!
//! The backend [`crate::curve::fiat::jubjub_scalar_64`] is produced by
//! fiat-crypto, exactly like the other `*_scalar_64` modules. Regenerate it with:
//!
//! ```text
//! word_by_word_montgomery --lang Rust --inline jubjub_scalar 64 \
//!   0x0e7db4ea6533afa906673b0101343b00a6682093ccc81082d0970e5ed6f72cb7 \
//!   mul square add sub opp from_montgomery to_montgomery nonzero selectznz \
//!   to_bytes from_bytes one msat divstep divstep_precomp
//! ```

use crate::curve::fiat::jubjub_scalar_64::*;
use crate::curve::field::{Field, Sign};
use crate::fiat_field_montgomery_impl;
use crate::mp::ct::{Choice, CtEqual, CtZero};
use crate::params::jubjub::ORDER_LIMBS;

/// Number of 64-bit limbs of a scalar field element (252 bits -> 4 limbs).
const GM_LIMBS_SIZE: usize = 4;

fiat_field_montgomery_impl!(
    #[doc = "Element of the Jubjub scalar prime field where r = 0x0e7db4ea6533afa906673b0101343b00a6682093ccc81082d0970e5ed6f72cb7"]
    Scalar,
    252,
    ORDER_LIMBS,
    GM_LIMBS_SIZE,
    fiat_jubjub_scalar_non_montgomery_domain_field_element,
    fiat_jubjub_scalar_nonzero,
    fiat_jubjub_scalar_add,
    fiat_jubjub_scalar_sub,
    fiat_jubjub_scalar_mul,
    fiat_jubjub_scalar_square,
    fiat_jubjub_scalar_opp,
    fiat_jubjub_scalar_to_bytes,
    fiat_jubjub_scalar_from_bytes,
    fiat_jubjub_scalar_montgomery_domain_field_element,
    fiat_jubjub_scalar_to_montgomery,
    fiat_jubjub_scalar_from_montgomery,
    fiat_jubjub_scalar_selectznz,
    fiat_jubjub_scalar_msat,
    fiat_jubjub_scalar_divstep,
    fiat_jubjub_scalar_divstep_precomp,
    le
);

impl Scalar {
    /// Get the multiplicative inverse.
    ///
    /// This currently delegates to the representation-agnostic Bernstein-Yang
    /// "safegcd" inversion ([`Self::inverse_safegcd`]); see the note on
    /// [`crate::curve::bls12_381::fp::Fp::inverse`].
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

    #[test]
    fn modulus_is_r() {
        // r - 1, big-endian
        const R_MINUS_1_BE: [u8; 32] = [
            0x0e, 0x7d, 0xb4, 0xea, 0x65, 0x33, 0xaf, 0xa9, 0x06, 0x67, 0x3b, 0x01, 0x01, 0x34,
            0x3b, 0x00, 0xa6, 0x68, 0x20, 0x93, 0xcc, 0xc8, 0x10, 0x82, 0xd0, 0x97, 0x0e, 0x5e,
            0xd6, 0xf7, 0x2c, 0xb6,
        ];
        let minus_one = Scalar::zero() - Scalar::one();
        assert_eq!(minus_one.to_bytes_be(), R_MINUS_1_BE);
    }
}
