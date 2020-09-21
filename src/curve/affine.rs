use super::weierstrass::WeierstrassCurve;
use crate::curve::field::{Field, FieldSqrt, Sign};
use core::ops::{Add, Mul, Sub};

/// Affine point operation over Field element FE
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Point<FE> {
    pub x: FE,
    pub y: FE,
}

impl<FE: Field> Point<FE> {
    pub fn to_coordinate(&self) -> (&FE, &FE) {
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
    pub fn decompress<C: WeierstrassCurve<FieldElement = FE>>(
        x: &FE,
        y_sign: Sign,
        curve: C,
    ) -> Option<Self> {
        // Y^2 = X^3 - A*X + b
        let yy = x.square() * x + (curve.a() * &x) + curve.b();
        match yy.sqrt().into_option() {
            None => None,
            Some(y) => {
                let found_sign = y.sign();
                let ny = -y.clone();
                if found_sign == y_sign {
                    Some(Point { x: x.clone(), y })
                } else {
                    Some(Point {
                        x: x.clone(),
                        y: ny,
                    })
                }
            }
        }
    }
}

impl<FE> Point<FE>
where
    FE: Field,
    // extend field operation to `&FE OP &FE`
    for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
    for<'a, 'b> &'a FE: Sub<&'b FE, Output = FE>,
{
    pub fn from_coordinate<C: WeierstrassCurve<FieldElement = FE>>(
        x: &FE,
        y: &FE,
        curve: C,
    ) -> Option<Self> {
        let y2 = y.square();
        let x3 = x.square() * x;
        let ax = curve.a() * x;

        if y2 == x3 + ax + curve.b() {
            Some(Point {
                x: x.clone(),
                y: y.clone(),
            })
        } else {
            None
        }
    }

    pub fn double<C: WeierstrassCurve<FieldElement = FE>>(&self, curve: C) -> Self {
        let Point {
            x: ref x1,
            y: ref y1,
        } = self;
        let l = (FE::from(3u64) * (x1.square()) + curve.a()) * (y1.double()).inverse();
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
