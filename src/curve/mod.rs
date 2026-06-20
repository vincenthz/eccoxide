//! Generic Curve mathematic and some known specific curve implementation
//!
//! For the generic mathematic:
//! * field: Field abstractions
//! * affine: Affine point on short weierstrass curve
//! * projective: Projective point on short weierstrass curve
//! * weierstrass: Abstraction for short weierstrass curve
//! * montgomery: Abstraction for montgomery curve
//! * edwards: Abstraction for (twisted) edwards curve
//!
//! For implementation of specific curve:
//! * sec2 (e.g. p192r1, p5p256k1, p256k1, p384r1, p521r1)

#[cfg(any(
    feature = "p112r1",
    feature = "p112r2",
    feature = "p128r1",
    feature = "p128r2",
    feature = "p160k1",
    feature = "p160r1",
    feature = "p160r2",
))]
pub(crate) mod bigint; // module used for compat and naive implementations

pub(crate) mod fiat;

pub mod affine;
pub mod edwards;
pub mod field;
pub mod montgomery;
pub mod projective;
pub mod weierstrass;

pub use field::Sign;

// exports the SEC2 curves
pub mod sec2;

// curve25519 / edwards25519
#[cfg(feature = "curve25519")]
pub mod curve25519;

// curve448 ("Goldilocks")
#[cfg(feature = "curve448")]
pub mod curve448;
