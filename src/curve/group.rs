//! Group abstraction for elliptic curve points
//!
//! A [`CurveGroup`] is a prime-order group written additively, whose elements
//! are curve points: it has an identity element (the point at infinity for
//! weierstrass curves, the neutral element otherwise), a generator, point
//! addition/negation, and scalar multiplication by elements of the associated
//! [`Field`] of scalars (the prime field of order the group order).
//!
//! This is implemented by the concrete point types of each enabled curve
//! (e.g. `curve::sec2::p256r1::Point`, `curve::curve25519::Point`,
//! `curve::curve25519::ristretto255::RistrettoPoint`), allowing protocols to
//! be written generically over any curve group.

use crate::curve::field::Field;
use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

/// Abstract trait for a prime-order group of points on an elliptic curve
///
/// The group law is written additively: `+` adds two points (handling
/// doubling and the identity), `-` negates or subtracts, and constant-time
/// scalar multiplication is the `*` operator with a [`CurveGroup::Scalar`]
/// reference.
///
/// All the provided operations are expected to be constant-time with respect
/// to secret data, except [`CurveGroup::mul_vartime`] which explicitly trades
/// that property for speed and must only be used with public scalars.
pub trait CurveGroup:
    Sized
    + 'static
    + Send
    + Sync
    + Clone
    + PartialEq
    + Eq
    + fmt::Debug
    + Neg<Output = Self>
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a <Self as CurveGroup>::Scalar, Output = Self>
{
    /// The field of scalars acting on the group: the prime field of order the
    /// group order `n`, so that for any point `P`, `n * P = IDENTITY`.
    type Scalar: Field;

    /// The neutral element of the group which is the additive zero: `P + IDENTITY = P`
    const IDENTITY: Self;

    /// The group generator (base point)
    const GENERATOR: Self;

    /// Double the point: `2 * self`, equivalent to `self + self` but usually
    /// implemented with a faster dedicated formula
    fn double(&self) -> Self;

    /// Constant-time fixed-base scalar multiplication: `scalar * GENERATOR`.
    ///
    /// Implementations typically use a precomputed table for the generator,
    /// making this substantially faster than `&Self::GENERATOR * scalar`.
    fn mul_base(scalar: &Self::Scalar) -> Self;

    /// Variable-time scalar multiplication: `scalar * self`.
    ///
    /// Faster than the constant-time `*` operator where implemented, but its
    /// running time depends on the scalar value; only use it when the scalar
    /// is public.
    fn mul_vartime(&self, scalar: &Self::Scalar) -> Self;
}

#[cfg(test)]
mod tests {
    use super::CurveGroup;

    /// Check the group laws and the consistency of the trait operations for
    /// one implementation
    fn group_laws<G: CurveGroup>() {
        let g = G::GENERATOR;
        let id = G::IDENTITY;

        // the identity is the additive zero
        assert_eq!(g.clone() + &id, g);
        assert_eq!(id.clone() + &g, g);

        // negation and subtraction return to the identity
        assert_eq!(g.clone() - &g, id);
        assert_eq!(g.clone() + (-g.clone()), id);

        // doubling is addition with self
        let g2 = g.double();
        assert_eq!(g2, g.clone() + &g);

        // the scalar multiplications match repeated addition
        let three = G::Scalar::from(3);
        let g3 = g2.clone() + &g;
        assert_eq!(g.clone() * &three, g3);
        assert_eq!(G::mul_base(&three), g3);
        assert_eq!(g.mul_vartime(&three), g3);
        assert_eq!(g.clone() * &three, g3);

        // scalar multiplication is linear: (a + b) * P == a*P + b*P
        let a = G::Scalar::from(0x1234_5678_9abc_def0);
        let b = G::Scalar::from(0x0fed_cba9_8765_4321);
        let ab = a.clone() + &b;
        assert_eq!(G::mul_base(&ab), G::mul_base(&a) + G::mul_base(&b));
        assert_eq!(g.clone() * &ab, g.clone() * &a + g.clone() * &b);
    }

    // one a=-3 curve, one a=0 curve, and the edwards / ristretto groups
    #[cfg(feature = "p256r1")]
    #[test]
    fn p256r1() {
        group_laws::<crate::curve::sec2::p256r1::Point>()
    }

    #[cfg(feature = "p256k1")]
    #[test]
    fn p256k1() {
        group_laws::<crate::curve::sec2::p256k1::Point>()
    }

    #[cfg(feature = "curve25519")]
    #[test]
    fn edwards25519() {
        group_laws::<crate::curve::curve25519::Point>()
    }

    #[cfg(feature = "bls12-381")]
    #[test]
    fn bls12_381_g1() {
        group_laws::<crate::curve::bls12_381::g1::Point>()
    }

    #[cfg(feature = "bls12-381")]
    #[test]
    fn bls12_381_g2() {
        group_laws::<crate::curve::bls12_381::g2::Point>()
    }

    #[cfg(feature = "jubjub")]
    #[test]
    fn jubjub() {
        group_laws::<crate::curve::jubjub::Point>()
    }

    #[cfg(feature = "ristretto255")]
    #[test]
    fn ristretto255() {
        group_laws::<crate::curve::curve25519::ristretto255::RistrettoPoint>()
    }
}
