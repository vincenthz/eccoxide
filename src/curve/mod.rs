//! Generic Curve mathematic and some known specific curve implementation
//!
//! For the generic mathematic:
//! * field: Field abstractions
//! * affine: Affine point on short weierstrass curve
//! * projective: Projective point on short weierstrass curve
//! * weierstrass: Abstraction for short weierstrass curve
//!
//! For implementation of specific curve:
//! * sec2 (e.g. p192r1, p5p256k1, p256k1, p384r1, p521r1)

pub(crate) mod bigint; // module used for compat and naive implementations
pub(crate) mod fiat;

pub mod affine;
pub mod field;
pub mod projective;
pub mod weierstrass;

pub use field::Sign;

// exports the SEC2 curves
pub mod sec2;
