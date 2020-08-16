use crate::params::sec2::p256k1::*;

mod fe;
mod gm;

pub use fe::FieldElement;
pub use gm::Scalar;

/// Affine Point on the curve
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PointAffine {
    x: FieldElement,
    y: FieldElement,
}

/// Point on the curve
#[derive(Clone, Debug)]
pub struct Point {
    x: FieldElement,
    y: FieldElement,
    z: FieldElement,
}

impl PartialEq for Point {
    fn eq(&self, other: &Point) -> bool {
        self.to_affine() == other.to_affine()
    }
}

impl Eq for Point {}

lazy_static! {
    static ref A: FieldElement = FieldElement::from_bytes(&A_BYTES).unwrap();
    static ref B: FieldElement = FieldElement::from_bytes(&B_BYTES).unwrap();
    static ref B3: FieldElement = FieldElement::from_bytes(&B3_BYTES).unwrap();
    static ref GX: FieldElement = FieldElement::from_bytes(&GX_BYTES).unwrap();
    static ref GY: FieldElement = FieldElement::from_bytes(&GY_BYTES).unwrap();
}

impl PointAffine {
    /// Curve generator point
    pub fn generator() -> Self {
        PointAffine {
            x: GX.clone(),
            y: GY.clone(),
        }
    }

    // check if y^2 = x^3 + a*x + b (mod p) holds
    pub fn from_coordinate(x: &FieldElement, y: &FieldElement) -> Option<Self> {
        let y2 = y * y;
        let x3 = x.power(3);
        let ax = &*A * x;

        if y2 == x3 + ax + &*B {
            Some(PointAffine {
                x: x.clone(),
                y: y.clone(),
            })
        } else {
            None
        }
    }

    pub fn to_coordinate(&self) -> (&FieldElement, &FieldElement) {
        (&self.x, &self.y)
    }

    pub fn double(&self) -> PointAffine {
        let PointAffine {
            x: ref x1,
            y: ref y1,
        } = self;
        let l = (FieldElement::from_u64(3) * (x1 * x1) + &*A)
            * (FieldElement::from_u64(2) * y1).inverse();
        let l2 = &l * &l;
        let x3 = l2 - FieldElement::from_u64(2) * x1;
        let y3 = l * (x1 - &x3) - y1;
        PointAffine { x: x3, y: y3 }
    }

    pub fn compress(&self) -> (&FieldElement, bool) {
        (&self.x, self.y.high_bit_set())
    }

    pub fn decompress(x: &FieldElement, bit: bool) -> Option<Self> {
        // Y^2 = X^3 - A*X + b
        let x = x.clone();
        //let x2 = x.square();
        //let x3 = &x * &x2;
        let yy = x.power(3) + (&*A * &x) + &*B;
        match yy.sqrt().into_option() {
            None => None,
            Some(y) => {
                let bit_set = y.high_bit_set();
                let ny = -&y;
                if bit == bit_set {
                    Some(PointAffine { x, y })
                } else {
                    Some(PointAffine { x, y: ny })
                }
            }
        }
    }
}

impl<'a, 'b> std::ops::Add<&'b PointAffine> for &'a PointAffine {
    type Output = PointAffine;
    fn add(self, other: &'b PointAffine) -> PointAffine {
        let PointAffine {
            x: ref x1,
            y: ref y1,
        } = &self;
        let PointAffine {
            x: ref x2,
            y: ref y2,
        } = &other;
        let l = (y1 - y2) * (x1 - x2).inverse();
        let l2 = &l * &l;
        let x3 = l2 - x1 - x2;
        let y3 = &l * (x1 - &x3) - y1;
        PointAffine { x: x3, y: y3 }
    }
}

impl Point {
    /// Curve generator point
    pub fn generator() -> Self {
        Point {
            x: GX.clone(),
            y: GY.clone(),
            z: FieldElement::one(),
        }
    }

    /// Point at infinity
    pub fn infinity() -> Self {
        Point {
            x: FieldElement::zero(),
            y: FieldElement::one(),
            z: FieldElement::zero(),
        }
    }

    pub fn from_affine(p: &PointAffine) -> Self {
        Point {
            x: p.x.clone(),
            y: p.y.clone(),
            z: FieldElement::one(),
        }
    }

    pub fn to_affine(&self) -> Option<PointAffine> {
        if self.z == FieldElement::one() {
            return Some(PointAffine {
                x: self.x.clone(),
                y: self.y.clone(),
            });
        }
        if self.z.is_zero() {
            None
        } else {
            let inv = self.z.inverse();
            Some(PointAffine {
                x: &self.x * &inv,
                y: &self.y * &inv,
            })
        }
    }

    pub fn normalize(&mut self) {
        let zinv = self.z.inverse();

        self.x = &self.x * &zinv;
        self.y = &self.y * &zinv;
        self.z = FieldElement::one()
    }

    fn add_different<'b>(&self, other: &'b Point) -> Point {
        assert!(self != other);

        let Point {
            x: ref x1,
            y: ref y1,
            z: ref z1,
        } = &self;
        let Point {
            x: ref x2,
            y: ref y2,
            z: ref z2,
        } = &other;

        // Algorithm 1 from (1) - addition formula for arbitrary a
        //
        // ```magma
        // ADD := function ( X1 , Y1 , Z1 , X2 , Y2 , Z2 ,a , b3 )
        //     t0 := X1 * X2 ; t1 := Y1 * Y2 ; t2 := Z1 * Z2 ;
        //     t3 := X1 + Y1 ; t4 := X2 + Y2 ; t3 := t3 * t4 ;
        //     t4 := t0 + t1 ; t3 := t3 - t4 ; t4 := X1 + Z1 ;
        //     t5 := X2 + Z2 ; t4 := t4 * t5 ; t5 := t0 + t2 ;
        //     t4 := t4 - t5 ; t5 := Y1 + Z1 ; X3 := Y2 + Z2 ;
        //     t5 := t5 * X3 ; X3 := t1 + t2 ; t5 := t5 - X3 ;
        //     Z3 := a * t4 ; X3 := b3 * t2 ; Z3 := X3 + Z3 ;
        //     X3 := t1 - Z3 ; Z3 := t1 + Z3 ; Y3 := X3 * Z3 ;
        //     t1 := t0 + t0 ; t1 := t1 + t0 ; t2 := a * t2 ;
        //     t4 := b3 * t4 ; t1 := t1 + t2 ; t2 := t0 - t2 ;
        //     t2 := a * t2 ; t4 := t4 + t2 ; t0 := t1 * t4 ;
        //     Y3 := Y3 + t0 ; t0 := t5 * t4 ; X3 := t3 * X3 ;
        //     X3 := X3 - t0 ; t0 := t3 * t1 ; Z3 := t5 * Z3 ;
        //     Z3 := Z3 + t0 ;
        //     return X3 , Y3 , Z3 ;
        // end function ;
        // ```

        let t0 = x1 * x2;
        let t1 = y1 * y2;
        let t2 = z1 * z2;
        let t3 = x1 + y1;
        let t4 = x2 + y2;
        let t3 = t3 * t4;
        let t4 = &t0 + &t1;
        let t3 = t3 - &t4;
        let t4 = x1 + z1;
        let t5 = x2 + z2;
        let t4 = t4 * &t5;
        let t5 = &t0 + &t2;
        let t4 = t4 - &t5;
        let t5 = y1 + z1;
        let x3 = y2 + z2;
        let t5 = t5 * &x3;
        let x3 = &t1 + &t2;
        let t5 = t5 - &x3;
        let z3 = &*A * &t4;
        let x3 = &*B3 * &t2;
        let z3 = &x3 + z3;
        let x3 = &t1 - &z3;
        let z3 = &t1 + &z3;
        let y3 = &x3 * &z3;
        let t1 = t0.double();
        let t1 = t1 + &t0;
        let t2 = &*A * t2;
        let t4 = &*B3 * &t4;
        let t1 = t1 + &t2;
        let t2 = &t0 - t2;
        let t2 = &*A * t2;
        let t4 = t4 + &t2;
        let t0 = &t1 * &t4;
        let y3 = y3 + t0;
        let t0 = &t5 * &t4;
        let x3 = &t3 * x3;
        let x3 = x3 - &t0;
        let t0 = &t3 * &t1;
        let z3 = &t5 * &z3;
        let z3 = z3 + t0;

        Point {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    pub fn double(&self) -> Self {
        let Point {
            ref x,
            ref y,
            ref z,
        } = &self;

        // Algorithm 3 from (1) - doubling formula for arbitrary a
        // ```magma
        // DBL := function (X ,Y ,Z ,a , b3 )
        //    t0 := X ^2; t1 := Y ^2; t2 := Z ^2;
        //    t3 := X * Y ; t3 := t3 + t3 ; Z3 := X * Z ;
        //    Z3 := Z3 + Z3 ; X3 := a * Z3 ; Y3 := b3 * t2 ;
        //    Y3 := X3 + Y3 ; X3 := t1 - Y3 ; Y3 := t1 + Y3 ;
        //    Y3 := X3 * Y3 ; X3 := t3 * X3 ; Z3 := b3 * Z3 ;
        //    t2 := a * t2 ; t3 := t0 - t2 ; t3 := a * t3 ;
        //    t3 := t3 + Z3 ; Z3 := t0 + t0 ; t0 := Z3 + t0 ;
        //    t0 := t0 + t2 ; t0 := t0 * t3 ; Y3 := Y3 + t0 ;
        //    t2 := Y * Z ; t2 := t2 + t2 ; t0 := t2 * t3 ;
        //    X3 := X3 - t0 ; Z3 := t2 * t1 ; Z3 := Z3 + Z3 ;
        //    Z3 := Z3 + Z3 ;
        //    return X3 , Y3 , Z3 ;
        // end function ;
        // ```
        let t0 = x * x;
        let t1 = y * y;
        let t2 = z * z;
        let t3 = x * y;
        let t3 = t3.double();
        let z3 = x * z;
        let z3 = &z3 + &z3;
        let x3 = &*A * &z3;
        let y3 = &*B3 * &t2;
        let y3 = &x3 + &y3;
        let x3 = &t1 - &y3;
        let y3 = &t1 + &y3;
        let y3 = &x3 * &y3;
        let x3 = &t3 * &x3;
        let z3 = &*B3 * &z3;
        let t2 = &*A * &t2;
        let t3 = &t0 - &t2;
        let t3 = &*A * &t3;
        let t3 = &t3 + &z3;
        let z3 = &t0 + &t0;
        let t0 = &z3 + &t0;
        let t0 = &t0 + &t2;
        let t0 = &t0 * &t3;
        let y3 = &y3 + &t0;
        let t2 = y * z;
        let t2 = &t2 + &t2;
        let t0 = &t2 * &t3;
        let x3 = &x3 - &t0;
        let z3 = &t2 * &t1;
        let z3 = &z3 + &z3;
        let z3 = &z3 + &z3;

        Point {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    fn add_or_double<'b>(&self, other: &'b Point) -> Point {
        if self == other {
            self.double()
        } else {
            self.add_different(other)
        }
    }

    /// scalar multiplication : `n * self` with double-and-add algorithm with increasing index
    fn scalar_mul_daa_limbs8(&self, n: &[u8]) -> Self {
        let mut a = self.clone();
        let mut q = Point::infinity();

        for digit in n.iter().rev() {
            for i in 0..8 {
                if digit & (1 << i) != 0 {
                    q = &q + &a;
                }
                a = a.double()
            }
        }
        q
    }

    /*
    fn xx() {
        let x1 = self;
        let x8 = x1.square_rep(8)
    }
    */
}

impl From<PointAffine> for Point {
    fn from(p: PointAffine) -> Self {
        Point {
            x: p.x,
            y: p.y,
            z: FieldElement::one(),
        }
    }
}

impl From<&PointAffine> for Point {
    fn from(p: &PointAffine) -> Self {
        Point::from_affine(p)
    }
}

// *************
// Point Negation
// *************

impl std::ops::Neg for Point {
    type Output = Point;

    fn neg(self) -> Self::Output {
        Point {
            x: self.x,
            y: -self.y,
            z: self.z,
        }
    }
}

impl<'a> std::ops::Neg for &'a Point {
    type Output = Point;

    fn neg(self) -> Self::Output {
        Point {
            x: self.x.clone(),
            y: -&self.y,
            z: self.z.clone(),
        }
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
        self.scalar_mul_daa_limbs8(&other.to_bytes())
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
