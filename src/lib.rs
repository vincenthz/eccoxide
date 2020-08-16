//! Elliptic Curve Cryptography

#[macro_use]
extern crate lazy_static;

pub mod curve;
pub(crate) mod mp;
pub mod params;

#[cfg(test)]
mod tests;
