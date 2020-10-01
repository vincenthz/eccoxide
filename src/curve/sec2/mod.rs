//! Prime Elliptic Curve defined in [SEC2](https://www.secg.org/SEC2-Ver-1.0.pdf)
mod bigint;

#[cfg(feature = "p256k1")]
pub mod p256k1;
#[cfg(feature = "p256r1")]
pub mod p256r1;
#[cfg(feature = "p384r1")]
pub mod p384r1;
#[cfg(feature = "p521r1")]
pub mod p521r1;

pub use bigint::*;
