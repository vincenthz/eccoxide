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
//! let generator = Point::generator();
//! let same_generator = &generator + Point::infinity();
//!
//! // transform the point to affine coordinate
//! let point = generator.to_affine().unwrap();
//! let (x, y) = point.to_coordinate();
//! ```

#[macro_use]
extern crate lazy_static;

pub mod curve;
pub(crate) mod mp;
pub mod params;

#[cfg(test)]
mod tests;
