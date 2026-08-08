//! BLS12-381
//!
//! A pairing-friendly curve widely used for BLS signatures and zk-SNARKs. It is
//! built around two prime fields:
//!
//! * the base field [`Fp`](fp::Fp) of order `p` (381-bit), and
//! * the scalar field [`Scalar`](scalar::Scalar) of order `r` (255-bit), which
//!   is the order of the prime-order subgroups.
//!
//! On top of those, [`g1`] provides the G1 point group over `Fp` and [`g2`] the
//! G2 point group over the quadratic extension [`Fp2`](fp2::Fp2), both with a
//! [`CurveGroup`](crate::curve::group::CurveGroup) implementation.
//!
//! The full extension tower ([`Fp2`](fp2::Fp2), [`Fp6`](fp6::Fp6),
//! [`Fp12`](fp12::Fp12)) is built up on top of the base field, and [`pairing`]
//! provides the optimal-ate pairing `e: G1 × G2 → Fp12`, both as the single
//! [`pairing`](pairing::pairing) function and split into its
//! [`miller_loop`](pairing::miller_loop) /
//! [`final_exponentiation`](pairing::MillerLoopResult::final_exponentiation)
//! halves, so that a product of pairings costs one final exponentiation.
//!
//! Points of either group serialize to the standard zcash/IETF compressed and
//! uncompressed encodings through the `to_compressed` / `from_compressed` /
//! `to_uncompressed` / `from_uncompressed` methods of their `Point` and
//! `PointAffine` types. Decoding validates membership of the prime-order
//! subgroup, which the curve equation alone does not imply; the
//! `_oncurve_only` decoders and the `is_in_subgroup` predicate expose the two
//! halves of that separately.

pub mod fp;
pub mod fp12;
pub mod fp2;
pub mod fp6;
pub mod g1;
pub mod g2;
#[cfg(feature = "bls12-381-hash-to-curve")]
pub mod hash_to_curve;
pub mod pairing;
pub mod scalar;
mod serialize;

/// Absolute value of the BLS12-381 seed parameter `x = -0xd201000000010000`,
/// from which `p` and `r` are derived (`r = x⁴ - x² + 1`).
///
/// The seed drives the Miller loop of the [`pairing`] and the exponentiations of
/// its final step, and — since `p ≡ x (mod r)` — the subgroup membership checks
/// of [`g1`] and [`g2`] as well. Its sign is carried separately in
/// [`BLS_X_IS_NEGATIVE`] rather than in the type, since every use of it iterates
/// over the bits of the magnitude.
pub(crate) const BLS_X: u64 = 0xd201000000010000;

/// Whether the seed [`BLS_X`] stands for a negative parameter (it does).
pub(crate) const BLS_X_IS_NEGATIVE: bool = true;

pub use fp::Fp;
pub use fp12::Fp12;
pub use fp2::Fp2;
pub use fp6::Fp6;
pub use scalar::Scalar;
