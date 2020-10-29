//! Prime Elliptic Curve defined in [SEC2](https://www.secg.org/sec2-v2.pdf)

#[cfg(any(
    feature = "p112r1",
    feature = "p112r2",
    feature = "p128r1",
    feature = "p128r2",
    feature = "p160k1",
    feature = "p160r1",
    feature = "p160r2",
))]
mod bigint;

#[cfg(feature = "p192k1")]
pub mod p192k1;
#[cfg(feature = "p192r1")]
pub mod p192r1;
#[cfg(feature = "p224k1")]
pub mod p224k1;
#[cfg(feature = "p224r1")]
pub mod p224r1;
#[cfg(feature = "p256k1")]
pub mod p256k1;
#[cfg(feature = "p256r1")]
pub mod p256r1;
#[cfg(feature = "p384r1")]
pub mod p384r1;
#[cfg(feature = "p521r1")]
pub mod p521r1;

#[cfg(any(
    feature = "p112r1",
    feature = "p112r2",
    feature = "p128r1",
    feature = "p128r2",
    feature = "p160k1",
    feature = "p160r1",
    feature = "p160r2",
))]
pub use bigint::*;
