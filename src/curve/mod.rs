pub(crate) mod bigint; // module used for compat and naive implementations
pub(crate) mod fiat;

pub mod affine;
pub mod field;
pub mod projective;
pub mod weierstrass;

// exports the SEC2 curves
pub mod sec2;
