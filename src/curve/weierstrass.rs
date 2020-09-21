//! short ℘ curves

/// Weierstrass curve are defined as y^{2} = x^{3} + Ax + B
pub trait WeierstrassCurve: Copy + Clone {
    type FieldElement;

    // Weirstrass A parameter
    fn a(self) -> &'static Self::FieldElement;
    // Weirstrass B parameter
    fn b(self) -> &'static Self::FieldElement;
    // Weirstrsass B parameter multiplied by 3
    fn b3(self) -> &'static Self::FieldElement;
}

/// Weierstrass curves with with A=0
pub trait WeierstrassCurveA0: WeierstrassCurve {}

/// Weierstrass curves with with A=-3
pub trait WeierstrassCurveAM3: WeierstrassCurve {}
