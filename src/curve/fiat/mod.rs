#![allow(dead_code)]
#![allow(non_camel_case_types)]

#[cfg(feature = "bls12-381")]
pub mod bls12_381_64;
#[cfg(feature = "bls12-381")]
pub mod bls12_381_scalar_64;
#[cfg(feature = "curve25519")]
pub mod curve25519_64;
#[cfg(feature = "curve25519")]
pub mod curve25519_scalar_64;
#[cfg(feature = "jubjub")]
pub mod jubjub_scalar_64;
#[cfg(feature = "p192k1")]
pub mod p192k1_64;
#[cfg(feature = "p192k1")]
pub mod p192k1_scalar_64;
#[cfg(feature = "p192r1")]
pub mod p192r1_64;
#[cfg(feature = "p192r1")]
pub mod p192r1_scalar_64;
#[cfg(feature = "p224r1")]
pub mod p224_64;
#[cfg(feature = "p224r1")]
pub mod p224_scalar_64;
#[cfg(feature = "p224k1")]
pub mod p224k1_64;
#[cfg(feature = "p224k1")]
pub mod p224k1_scalar_64;
#[cfg(feature = "p256r1")]
pub mod p256_64;
#[cfg(feature = "p256r1")]
pub mod p256_scalar_64;
#[cfg(feature = "p384r1")]
pub mod p384_64;
#[cfg(feature = "p384r1")]
pub mod p384_scalar_64;
#[cfg(feature = "curve448")]
pub mod p448_solinas_64;
#[cfg(feature = "p521r1")]
pub mod p521_64;
#[cfg(feature = "p521r1")]
pub mod p521_scalar_64;
#[cfg(feature = "p256k1")]
pub mod secp256k1_montgomery_64;
#[cfg(feature = "p256k1")]
pub mod secp256k1_montgomery_scalar_64;
#[cfg(feature = "sm2")]
pub mod sm2_64;
#[cfg(feature = "sm2")]
pub mod sm2_scalar_64;

mod curve_macros;
mod field_macros;
