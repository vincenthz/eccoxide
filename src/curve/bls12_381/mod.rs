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

pub mod fp;
pub mod fp12;
pub mod fp2;
pub mod fp6;
pub mod g1;
pub mod g2;
pub mod pairing;
pub mod scalar;

pub use fp::Fp;
pub use fp12::Fp12;
pub use fp2::Fp2;
pub use fp6::Fp6;
pub use scalar::Scalar;
