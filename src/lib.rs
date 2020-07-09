//! Elliptic Curve Cryptography

#[macro_use]
extern crate lazy_static;

pub mod curve;
pub mod params;

#[cfg(test)]
mod tests;
