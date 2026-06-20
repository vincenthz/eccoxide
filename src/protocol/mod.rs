//! Higher-level cryptographic protocols built on curves
//!
//! * [`x25519`]: X25519 Diffie-Hellman key agreement (RFC 7748) built on [`crate::curve::curve25519`]
//! * [`ed25519`]: Ed25519 digital signatures (RFC 8032) built on [`crate::curve::curve25519`]
//!

#[cfg(feature = "ed25519")]
pub mod ed25519;
#[cfg(feature = "x25519")]
pub mod x25519;
#[cfg(feature = "x448")]
pub mod x448;
