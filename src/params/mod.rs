//! Constant related to known elliptic curves

pub mod sec2;

#[cfg(feature = "bls12-381")]
pub mod bls12_381;

#[cfg(feature = "jubjub")]
pub mod jubjub;

#[cfg(all(feature = "curve25519", feature = "table"))]
pub mod curve25519;
