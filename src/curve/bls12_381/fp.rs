//! BLS12-381 base field `Fp`
//!
//! `Fp` is the prime field of order `p`, the 381-bit modulus over which the
//! curve (and later the extension fields `Fp2`, `Fp6`, `Fp12`) are defined.
//!
//! The arithmetic is provided by the fiat-crypto Montgomery backend; this
//! module only wires it into the generic [`Field`] abstraction, plus the
//! `p = 3 (mod 4)` square root needed for point (de)compression.

use crate::curve::fiat::bls12_381_64::*;
use crate::curve::field::{Field, FieldSqrt, Sign};
use crate::mp::ct::{Choice, CtEqual, CtOption, CtZero};
use crate::params::bls12_381::{P_LIMBS, P_PLUS1_DIV4_BYTES};
use crate::{fiat_field_montgomery_impl, fiat_field_sqrt_define};

/// Number of 64-bit limbs of a base field element (381 bits -> 6 limbs).
const FE_LIMBS_SIZE: usize = 6;

fiat_field_montgomery_impl!(
    #[doc = "Element of the BLS12-381 base prime field Fp where p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab"]
    Fp,
    381,
    P_LIMBS,
    FE_LIMBS_SIZE,
    fiat_bls12_381_non_montgomery_domain_field_element,
    fiat_bls12_381_nonzero,
    fiat_bls12_381_add,
    fiat_bls12_381_sub,
    fiat_bls12_381_mul,
    fiat_bls12_381_square,
    fiat_bls12_381_opp,
    fiat_bls12_381_to_bytes,
    fiat_bls12_381_from_bytes,
    fiat_bls12_381_montgomery_domain_field_element,
    fiat_bls12_381_to_montgomery,
    fiat_bls12_381_from_montgomery,
    fiat_bls12_381_selectznz,
    fiat_bls12_381_msat,
    fiat_bls12_381_divstep,
    fiat_bls12_381_divstep_precomp
);
fiat_field_sqrt_define!(Fp);

impl Fp {
    /// Get the multiplicative inverse.
    ///
    /// This currently delegates to the representation-agnostic Bernstein-Yang
    /// "safegcd" inversion ([`Self::inverse_safegcd`]) rather than a
    /// hand-crafted Fermat addition chain: it is constant-time and correct for
    /// the 381-bit modulus, at the cost of some speed. A dedicated addition
    /// chain can be introduced later as an optimisation.
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a
    /// panic.
    pub fn inverse(&self) -> Self {
        self.inverse_safegcd()
    }

    /// Compute the square root `x` of the field element such that `x*x = self`.
    ///
    /// Since `p = 3 (mod 4)`, the candidate root is `self^((p+1)/4)`; the
    /// returned [`CtOption`] is present only when that candidate actually
    /// squares back to `self` (i.e. `self` is a quadratic residue).
    pub fn sqrt(&self) -> CtOption<Self> {
        let candidate = self.power(&P_PLUS1_DIV4_BYTES);
        let check = candidate.square();
        CtOption::from((CtEqual::ct_eq(&check, self), candidate))
    }
}

#[cfg(test)]
mod tests {
    use super::Fp;
    use crate::{fiat_field_sqrt_unittest, fiat_field_unittest};

    fiat_field_unittest!(Fp);
    crate::fiat_field_safegcd_unittest!(Fp);
    fiat_field_sqrt_unittest!(Fp);
}
