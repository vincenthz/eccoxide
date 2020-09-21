use crate::curve::{
    affine,
    field::Sign,
    projective,
    weierstrass::{WeierstrassCurve, WeierstrassCurveA0},
};
use crate::params::sec2::p256k1::*;

mod fe;
mod gm;

pub use fe::FieldElement;
pub use gm::Scalar;

/// Affine Point on the curve
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PointAffine(affine::Point<FieldElement>);

/// Point on the curve
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Point(projective::Point<FieldElement>);

lazy_static! {
    static ref A: FieldElement = FieldElement::from_bytes(&A_BYTES).unwrap();
    static ref B: FieldElement = FieldElement::from_bytes(&B_BYTES).unwrap();
    static ref B3: FieldElement = FieldElement::from_bytes(&B3_BYTES).unwrap();
    static ref GX: FieldElement = FieldElement::from_bytes(&GX_BYTES).unwrap();
    static ref GY: FieldElement = FieldElement::from_bytes(&GY_BYTES).unwrap();
    static ref ORDER: FieldElement = FieldElement::from_bytes(&ORDER_BYTES).unwrap();
}

#[derive(Debug, Clone, Copy)]
pub struct Curve;

impl Curve {
    pub fn group_order(self) -> &'static FieldElement {
        &ORDER
    }
    pub fn generator() -> (&'static FieldElement, &'static FieldElement) {
        (&GX, &GY)
    }
}

impl WeierstrassCurve for Curve {
    type FieldElement = FieldElement;

    fn a(self) -> &'static Self::FieldElement {
        &A
    }

    fn b(self) -> &'static Self::FieldElement {
        &B
    }

    fn b3(self) -> &'static Self::FieldElement {
        &B3
    }
}

impl WeierstrassCurveA0 for Curve {}

impl PointAffine {
    /// Curve generator point
    pub fn generator() -> Self {
        PointAffine(affine::Point {
            x: GX.clone(),
            y: GY.clone(),
        })
    }

    // check if y^2 = x^3 + a*x + b (mod p) holds
    pub fn from_coordinate(x: &FieldElement, y: &FieldElement) -> Option<Self> {
        affine::Point::from_coordinate(x, y, Curve).map(PointAffine)
    }

    pub fn to_coordinate(&self) -> (&FieldElement, &FieldElement) {
        (&self.0.x, &self.0.y)
    }

    pub fn double(&self) -> PointAffine {
        PointAffine(affine::Point::double(&self.0, Curve))
    }

    pub fn compress(&self) -> (&FieldElement, Sign) {
        self.0.compress()
    }

    pub fn decompress(x: &FieldElement, sign: Sign) -> Option<Self> {
        affine::Point::decompress(x, sign, Curve).map(PointAffine)
    }
}

impl<'a, 'b> std::ops::Add<&'b PointAffine> for &'a PointAffine {
    type Output = PointAffine;
    fn add(self, other: &'b PointAffine) -> PointAffine {
        PointAffine(&self.0 + &other.0)
    }
}

impl Point {
    /// Curve generator point
    pub fn generator() -> Self {
        Point(projective::Point {
            x: GX.clone(),
            y: GY.clone(),
            z: FieldElement::one(),
        })
    }

    /// Point at infinity
    pub fn infinity() -> Self {
        Point(projective::Point::infinity())
    }

    pub fn from_affine(p: &PointAffine) -> Self {
        Point(projective::Point::from_affine(&p.0))
    }

    pub fn to_affine(&self) -> Option<PointAffine> {
        self.0.to_affine().map(PointAffine)
    }

    pub fn normalize(&mut self) {
        self.0.normalize()
    }

    fn add_or_double<'b>(&self, other: &'b Point) -> Point {
        Point(self.0.add_or_double_a0(&other.0, Curve))
    }
}

impl From<PointAffine> for Point {
    fn from(p: PointAffine) -> Self {
        Point(projective::Point::from_affine(&p.0))
    }
}

impl From<&PointAffine> for Point {
    fn from(p: &PointAffine) -> Self {
        Point(projective::Point::from_affine(&p.0))
    }
}

// *************
// Point Negation
// *************

impl std::ops::Neg for Point {
    type Output = Point;

    fn neg(self) -> Self::Output {
        Point(self.0.neg())
    }
}

impl<'a> std::ops::Neg for &'a Point {
    type Output = Point;

    fn neg(self) -> Self::Output {
        Point(self.0.clone().neg())
    }
}

// *************
// Point Scaling
// *************

// note that scalar multiplication is really defined for arbitrary scalar
// (of any size), not just the *field element* scalar defined in F(p).
// this semantic abuse makes it easier to use.

impl<'a, 'b> std::ops::Mul<&'b Scalar> for &'a Point {
    type Output = Point;

    fn mul(self, other: &'b Scalar) -> Point {
        Point(self.0.scale_a0(&other.to_bytes(), Curve))
    }
}

impl<'a, 'b> std::ops::Mul<&'b Point> for &'a Scalar {
    type Output = Point;

    fn mul(self, other: &'b Point) -> Point {
        other * self
    }
}

// **************
// Point Addition
// **************

impl<'a, 'b> std::ops::Add<&'b Point> for &'a Point {
    type Output = Point;

    fn add(self, other: &'b Point) -> Point {
        self.add_or_double(other)
    }
}

impl<'b> std::ops::Add<&'b Point> for Point {
    type Output = Point;

    fn add(self, other: &'b Point) -> Point {
        &self + other
    }
}

impl<'a> std::ops::Add<Point> for &'a Point {
    type Output = Point;

    fn add(self, other: Point) -> Point {
        self + &other
    }
}

impl std::ops::Add<Point> for Point {
    type Output = Point;

    fn add(self, other: Point) -> Point {
        &self + &other
    }
}

impl<'a, 'b> std::ops::Sub<&'b Point> for &'a Point {
    type Output = Point;

    fn sub(self, other: &'b Point) -> Point {
        self + (-other)
    }
}

impl std::ops::Sub<Point> for Point {
    type Output = Point;

    fn sub(self, other: Point) -> Point {
        &self - &other
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bytes() {
        let b: [u8; 32] = [
            1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 7, 8, 9, 0, 1,
            2, 3, 4,
        ];
        let s = Scalar::from_bytes(&b).unwrap();
        assert_eq!(b, s.to_bytes())
    }
    #[test]
    fn bytes_u64() {
        let b: [u8; 32] = [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0xce,
        ];

        let s1 = Scalar::from_bytes(&b).unwrap();
        let s2 = Scalar::from_u64(0xce);
        assert_eq!(s1, s2)
    }
}
