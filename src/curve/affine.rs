//! Affine Elliptic Curve Point defined over Field element as (X,Y)
//!
//! This only defined addition and multiplication operation, and coordinates decomposition
//!
//! Some other operations (negation, sub, etc) are also possible but this is not exhaustive
use super::weierstrass::WeierstrassCurve;
use crate::curve::field::{Field, FieldSqrt, Sign};
use crate::mp::ct::{Choice, CtOption};
use core::ops::{Add, Mul, Sub};

/// Affine point operation over Field element FE
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Point<FE> {
    pub x: FE,
    pub y: FE,
}

impl<FE: Field> Point<FE> {
    pub const fn to_coordinate(&self) -> (&FE, &FE) {
        (&self.x, &self.y)
    }

    pub fn compress(&self) -> (&FE, Sign) {
        (&self.x, self.y.sign())
    }
}

impl<FE> Point<FE>
where
    FE: FieldSqrt,
    // extend field operation to `&FE OP &FE`
    for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
    for<'a, 'b> &'a FE: Sub<&'b FE, Output = FE>,
{
    /// Recover the point from its x-coordinate and the requested sign of the
    /// y-coordinate (point decompression), in constant time.
    ///
    /// The recovery does not branch on the (secret) recovered value: the
    /// square-root presence is carried as a `Choice` rather than matched on,
    /// and the choice between `y` and `-y` is a branch-free constant-time
    /// select.
    ///
    /// The returned `CtOption` is present exactly when `x` is a valid
    /// compressed x-coordinate, i.e. when `x^3 + A*x + B` is a quadratic
    /// residue. When it is not present the carried point is a placeholder
    /// derived from the square-root candidate and must be discarded based on
    /// the presence choice.
    pub fn decompress<C: WeierstrassCurve<FieldElement = FE>>(
        x: &FE,
        y_sign: Sign,
    ) -> CtOption<Self> {
        let yy = x.square() * x + (&C::A * &x) + C::B;
        let (present, y) = yy.sqrt().into_parts();
        let ny = -y.clone();
        let matches = Choice::from(y.sign() == y_sign);
        let y = FE::ct_select(matches, &y, &ny);
        CtOption::from((present, Point { x: x.clone(), y }))
    }
}

impl<FE> Point<FE>
where
    FE: Field,
    // extend field operation to `&FE OP &FE`
    for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
    for<'a, 'b> &'a FE: Sub<&'b FE, Output = FE>,
{
    /// Create a point without verifying this is a valid point
    #[allow(unused)]
    pub(crate) const fn from_coordinate_unchecked<C: WeierstrassCurve<FieldElement = FE>>(
        x: FE,
        y: FE,
    ) -> Self {
        Point { x, y }
    }

    pub fn from_coordinate<C: WeierstrassCurve<FieldElement = FE>>(x: &FE, y: &FE) -> Option<Self> {
        let y2 = y.square();
        let x3 = x.square() * x;
        let ax = C::A * x;

        if y2 == x3 + ax + C::B {
            Some(Point {
                x: x.clone(),
                y: y.clone(),
            })
        } else {
            None
        }
    }

    pub fn double<C: WeierstrassCurve<FieldElement = FE>>(&self) -> Self {
        let Point {
            x: ref x1,
            y: ref y1,
        } = self;
        let l = (FE::from(3u64) * (x1.square()) + C::A) * (y1.double()).inverse();
        let l2 = l.square();
        let x3 = l2 - x1.double();
        let y3 = l * (x1 - &x3) - y1;
        Point { x: x3, y: y3 }
    }
}

impl<FE> Point<FE>
where
    FE: Field,
    // extend field operation to `&FE OP &FE`
    for<'a> &'a FE: Add<FE, Output = FE>,
    for<'a> &'a FE: Mul<FE, Output = FE>,
    for<'a> &'a FE: Sub<FE, Output = FE>,
    for<'a, 'b> &'a FE: Add<&'b FE, Output = FE>,
    for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
    for<'a, 'b> &'a FE: Sub<&'b FE, Output = FE>,
{
    pub fn add_different<'b>(&self, other: &'b Self) -> Self {
        let Point {
            x: ref x1,
            y: ref y1,
        } = &self;
        let Point {
            x: ref x2,
            y: ref y2,
        } = &other;
        let l = (y1 - y2) * (x1 - x2).inverse();
        let l2 = l.square();
        let x3 = l2 - x1 - x2;
        let y3 = l * (x1 - &x3) - y1;
        Point { x: x3, y: y3 }
    }
}

impl<'x, 'y, FE> std::ops::Add<&'y Point<FE>> for &'x Point<FE>
where
    FE: Field,
    for<'a> &'a FE: Add<FE, Output = FE>,
    for<'a> &'a FE: Mul<FE, Output = FE>,
    for<'a> &'a FE: Sub<FE, Output = FE>,
    for<'a, 'b> &'a FE: Add<&'b FE, Output = FE>,
    for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
    for<'a, 'b> &'a FE: Sub<&'b FE, Output = FE>,
{
    type Output = Point<FE>;
    fn add(self, other: &'y Point<FE>) -> Point<FE> {
        self.add_different(other)
    }
}
