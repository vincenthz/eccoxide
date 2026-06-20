//! (Twisted) Edwards curves
//!
//! All twisted edwards curve are defined as Ax^{2} + y^{2} = 1 + Dx^{2}y^{2}
//! The original (untwisted) edwards curve is the special case A=1

/// Twisted edwards curve are defined as Ax^{2} + y^{2} = 1 + Dx^{2}y^{2}
pub trait EdwardsCurve: Copy + Clone {
    type FieldElement;

    // Twisted edwards A parameter
    const A: Self::FieldElement;
    // Twisted edwards D parameter
    const D: Self::FieldElement;
    // Twisted edwards D parameter multiplied by 2 , used in the extended (X:Y:Z:T) addition
    const D2: Self::FieldElement;
}

/// Edwards curves with A=1 (the original, untwisted edwards form)
pub trait EdwardsCurveA1: EdwardsCurve {}

/// Edwards curves with A=-1 , which admit the fastest complete addition
pub trait EdwardsCurveAM1: EdwardsCurve {}
