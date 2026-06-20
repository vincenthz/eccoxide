//! Constant related to known elliptic curves

pub mod sec2;

#[cfg(all(feature = "curve25519", feature = "table"))]
pub mod curve25519;
