//! Elliptic Curve Cryptography
//!
//! Provide arithmetic functions to deal with
//! elliptic curve mathematics, from the finite field
//! arithmetics abstractions to point representation on
//! curves.
//!
//! It also provide those operations on multiples
//! standard curves, and use when possible
//! constant-time operations.
//!
//! Currently most lowlevel arithmetics is provided by
//! the fiat-crypto project
//!
//! Example to do basic arithmetic operation with this crate:
//!
//! ```
//! // use p256r1 for this example, but other curve available in sec2 module
//! // or future other hierarchy
//! use eccoxide::curve::sec2::p256r1::{FieldElement, Point};
//!
//! // addition in the underlying Field
//! let two = FieldElement::one() + FieldElement::one();
//!
//! // Get the generator and add the point at infinity
//! let same_generator = &Point::GENERATOR + Point::INFINITY;
//!
//! // transform the point to affine coordinate
//! let point = Point::GENERATOR.to_affine().unwrap();
//! let (x, y) = point.to_coordinate();
//! ```
//!
//! Using
//!
//! ```
//! // use p521r1 for this example
//! use eccoxide::curve::sec2::p521r1::{FieldElement, Scalar, Point};
//! use eccoxide::curve::Sign;
//!
//! // this is just unsecure example, but here it could be loading some secret from disk
//! let bytes : [u8;66] = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
//!
//! let secret_key = Scalar::from_bytes(&bytes).unwrap();
//! let public_key = &Point::GENERATOR * &secret_key;
//!
//! // serialize the public key to a standard-ish compress format for p521r1
//! let public_affine = public_key.to_affine().unwrap();
//! let (x, ysign) = public_affine.compress();
//! let format_byte : u8 = match ysign {
//!     Sign::Positive => 0x2,
//!     Sign::Negative => 0x3,
//! };
//!
//! let mut public_key_bytes = Vec::new();
//! public_key_bytes.push(format_byte);
//! public_key_bytes.extend_from_slice(&x.to_bytes());
//! ```

pub mod curve;
pub(crate) mod mp;
pub mod params;
#[cfg(any(feature = "x25519", feature = "ed25519"))]
pub mod protocol;

#[cfg(test)]
mod tests;
