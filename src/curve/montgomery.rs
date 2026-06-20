//! Montgomery curves
//!
//! All montgomery curve are defined as By^{2} = x^{3} + Ax^{2} + x

/// Montgomery curve are defined as By^{2} = x^{3} + Ax^{2} + x
pub trait MontgomeryCurve: Copy + Clone {
    type FieldElement;

    // Montgomery A parameter
    const A: Self::FieldElement;
    // Montgomery B parameter
    const B: Self::FieldElement;
    // (A + 2) / 4 , used in the differential (x-only) ladder
    const A24: Self::FieldElement;
}

/// Montgomery curves with B=1
pub trait MontgomeryCurveB1: MontgomeryCurve {}
