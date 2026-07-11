//! Higher-level cryptographic protocols built on curves
//!
//! * [`x25519`]: X25519 Diffie-Hellman key agreement (RFC 7748) built on [`crate::curve::curve25519`]
//! * [`ed25519`]: Ed25519 digital signatures (RFC 8032) built on [`crate::curve::curve25519`]
//! * [`ecdsa`]: ECDSA digital signatures (SEC1) built on the [`crate::curve::sec2`] curves
//!

#[cfg(feature = "ecdsa")]
pub mod ecdsa;
#[cfg(feature = "ed25519")]
pub mod ed25519;
#[cfg(feature = "x25519")]
pub mod x25519;
#[cfg(feature = "x448")]
pub mod x448;
