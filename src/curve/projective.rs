//! Projective Elliptic Curve Point defined over Field element as (X,Y,Z)
//!
//! This module implements the point addition formulas defined in
//! [Complete addition formulas for prime order elliptic curves](https://eprint.iacr.org/2015/1060.pdf)
//!
//! For more complete reading:
//!
//! * [Complete addition formulas for prime order elliptic curves](https://eprint.iacr.org/2015/1060.pdf) (1)
//! * Handbook of Elliptic and Hyperelliptic Curve Cryptography - Chapter 13
//! * [NIST.SP.800-186](https://csrc.nist.gov/publications/detail/sp/800-186/draft) : Appendix D & E

use super::affine;
use super::field::Field;
use super::weierstrass::{WeierstrassCurve, WeierstrassCurveA0};
use crate::mp::ct::{Choice, CtEqual};
use std::convert::TryFrom;
use std::ops::{Add, Mul, Neg, Sub};

/// Projective point with field element FE
///
/// Affine point associated with (X,Y,Z) : (X/Z, Y/Z)
///
/// Note that 2 points are equal if they are in the same equivalence class,
/// which is determined with 4 FieldElement multiplications.
///
/// Example: (1,2,1) and (2,4,2) are equal
#[derive(Clone, Debug)]
pub struct Point<FE> {
    pub x: FE,
    pub y: FE,
    pub z: FE,
}

impl<FE: Field + CtEqual> PartialEq for Point<FE>
where
    for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
{
    fn eq(&self, other: &Point<FE>) -> bool {
        self.is_equivalent(other).is_true()
    }
}

impl<FE: Field + CtEqual> Eq for Point<FE> where for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE> {}

impl<'y, FE: Field + CtEqual> CtEqual for Point<FE>
where
    for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
{
    fn ct_eq(&self, other: &Point<FE>) -> Choice {
        self.is_equivalent(other)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct AffineAtInfinity;

impl<FE: Field> TryFrom<Point<FE>> for affine::Point<FE>
where
    for<'a> &'a FE: Add<FE, Output = FE>,
    for<'a> &'a FE: Mul<FE, Output = FE>,
    for<'a> &'a FE: Sub<FE, Output = FE>,
    for<'a, 'b> &'a FE: Add<&'b FE, Output = FE>,
    for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
    for<'a, 'b> &'a FE: Sub<&'b FE, Output = FE>,
{
    type Error = AffineAtInfinity;

    fn try_from(p: Point<FE>) -> Result<affine::Point<FE>, Self::Error> {
        p.to_affine().ok_or(AffineAtInfinity)
    }
}

impl<FE> Point<FE>
where
    FE: Field + CtEqual,
    for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
{
    /// Check if a point are in the same equivalent class
    fn is_equivalent<'y>(&self, other: &Point<FE>) -> Choice {
        let nx1 = &self.x * &other.z;
        let nx2 = &other.x * &self.z;
        let ny1 = &self.y * &other.z;
        let ny2 = &other.y * &self.z;
        nx1.ct_eq(&nx2) & ny1.ct_eq(&ny2)
    }

    /// Check if a point is at infinity
    pub fn is_infinity(&self) -> Choice {
        self.z.ct_eq(&FE::zero())
    }
}

impl<FE> Point<FE>
where
    FE: Field,
{
    /// Returns the point at infinity
    pub fn infinity() -> Self {
        Point {
            x: FE::zero(),
            y: FE::one(),
            z: FE::zero(),
        }
    }

    pub fn from_affine(p: &affine::Point<FE>) -> Self {
        Point {
            x: p.x.clone(),
            y: p.y.clone(),
            z: FE::one(),
        }
    }
}

impl<FE> Point<FE>
where
    FE: Field,
    for<'a, 'b> &'a FE: Add<&'b FE, Output = FE>,
    for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
    for<'a, 'b> &'a FE: Sub<&'b FE, Output = FE>,
{
    pub fn normalize(&mut self) {
        if !self.z.is_zero() {
            let zinv = self.z.inverse();

            self.x = &self.x * &zinv;
            self.y = &self.y * &zinv;
            self.z = FE::one()
        }
    }
}

impl<FE: Field> Point<FE> {
    pub fn add_different<'x, 'y, C: WeierstrassCurve<FieldElement = FE>>(
        &'x self,
        other: &'y Point<FE>,
        curve: C,
    ) -> Point<FE>
    where
        for<'a> &'a FE: Add<FE, Output = FE>,
        for<'a> &'a FE: Mul<FE, Output = FE>,
        for<'a> &'a FE: Sub<FE, Output = FE>,
        for<'a, 'b> &'a FE: Add<&'b FE, Output = FE>,
        for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
        for<'a, 'b> &'a FE: Sub<&'b FE, Output = FE>,
    {
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

        let t0 = &self.x * &other.x;
        let t1 = &self.y * &other.y;
        let t2 = &self.z * &other.z;
        let t3 = &self.x + &self.y;
        let t4 = &other.x + &other.y;
        let t3 = t3 * t4;
        let t4 = &t0 + &t1;
        let t3 = t3 - &t4;
        let t4 = &self.x + &self.z;
        let t5 = &other.x + &other.z;
        let t4 = t4 * &t5;
        let t5 = &t0 + &t2;
        let t4 = t4 - &t5;
        let t5 = &self.y + &self.z;
        let x3 = &other.y + &other.z;
        let t5 = t5 * &x3;
        let x3 = &t1 + &t2;
        let t5 = t5 - &x3;
        let z3 = curve.a() * &t4;
        let x3 = curve.b3() * &t2;
        let z3 = &x3 + &z3;
        let x3 = &t1 - &z3;
        let z3 = &t1 + &z3;
        let y3 = &x3 * &z3;
        let t1 = t0.double();
        let t1 = t1 + &t0;
        let t2 = curve.a() * &t2;
        let t4 = curve.b3() * &t4;
        let t1 = t1 + &t2;
        let t2 = &t0 - &t2;
        let t2 = curve.a() * &t2;
        let t4 = t4 + &t2;
        let t0 = &t1 * &t4;
        let y3 = y3 + t0;
        let t0 = &t5 * &t4;
        let x3 = &t3 * &x3;
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

    pub fn add_different_a0<'x, 'y, C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveA0>(
        &'x self,
        other: &'y Point<FE>,
        curve: C,
    ) -> Point<FE>
    where
        for<'a> &'a FE: Add<FE, Output = FE>,
        for<'a> &'a FE: Mul<FE, Output = FE>,
        for<'a> &'a FE: Sub<FE, Output = FE>,
        for<'a, 'b> &'a FE: Add<&'b FE, Output = FE>,
        for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
        for<'a, 'b> &'a FE: Sub<&'b FE, Output = FE>,
    {
        //
        // Algorithm from (1) - addition formula for a=0
        //
        // ```magma
        // ADD := function ( X1 , Y1 , Z1 , X2 , Y2 , Z2 , b3 )
        //     t0 := X1 * X2 ; t1 := Y1 * Y2 ; t2 := Z1 * Z2 ;
        //     t3 := X1 + Y1 ; t4 := X2 + Y2 ; t3 := t3 * t4 ;
        //     t4 := t0 + t1 ; t3 := t3 - t4 ; t4 := Y1 + Z1 ;
        //     X3 := Y2 + Z2 ; t4 := t4 * X3 ; X3 := t1 + t2 ;
        //     t4 := t4 - X3 ; X3 := X1 + Z1 ; Y3 := X2 + Z2 ;
        //     X3 := X3 * Y3 ; Y3 := t0 + t2 ; Y3 := X3 - Y3 ;
        //     X3 := t0 + t0 ; t0 := X3 + t0 ; t2 := b3 * t2 ;
        //     Z3 := t1 + t2 ; t1 := t1 - t2 ; Y3 := b3 * Y3 ;
        //     X3 := t4 * Y3 ; t2 := t3 * t1 ; X3 := t2 - X3 ;
        //     Y3 := Y3 * t0 ; t1 := t1 * Z3 ; Y3 := t1 + Y3 ;
        //     t0 := t0 * t3 ; Z3 := Z3 * t4 ; Z3 := Z3 + t0 ;
        //     return X3 , Y3 , Z3 ;
        // end function ;
        // ```
        let t0 = &self.x * &other.x;
        let t1 = &self.y * &other.y;
        let t2 = &self.z * &other.z;
        let t3 = &self.x + &self.y;
        let t4 = &other.x + &other.y;
        let t3 = t3 * t4;
        let t4 = &t0 + &t1;
        let t3 = t3 - t4;
        let t4 = &self.y + &self.z;
        let x3 = &other.y + &other.z;
        let t4 = &t4 * &x3;
        let x3 = &t1 + &t2;
        let t4 = &t4 - &x3;
        let x3 = &self.x + &self.z;
        let y3 = &other.x + &other.z;
        let x3 = x3 * y3;
        let y3 = &t0 + &t2;
        let y3 = x3 - y3;
        let x3 = &t0 + &t0;
        let t0 = &x3 + &t0;
        let t2 = curve.b3() * &t2;
        let z3 = &t1 + &t2;
        let t1 = t1 - &t2;
        let y3 = curve.b3() * y3;
        let x3 = &t4 * &y3;
        let t2 = &t3 * &t1;
        let x3 = &t2 - x3;
        let y3 = y3 * &t0;
        let t1 = &t1 * &z3;
        let y3 = &t1 + y3;
        let t0 = t0 * t3;
        let z3 = z3 * t4;
        let z3 = z3 + t0;

        Point {
            x: x3,
            y: y3,
            z: z3,
        }
    }
}

impl<FE> Point<FE>
where
    FE: Field,
    for<'a> &'a FE: Add<FE, Output = FE>,
    for<'a> &'a FE: Mul<FE, Output = FE>,
    for<'a> &'a FE: Sub<FE, Output = FE>,
    for<'a, 'b> &'a FE: Add<&'b FE, Output = FE>,
    for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
    for<'a, 'b> &'a FE: Sub<&'b FE, Output = FE>,
{
    #[inline]
    pub fn double<C: WeierstrassCurve<FieldElement = FE>>(&self, curve: C) -> Self {
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
        let t0 = self.x.square();
        let t1 = self.y.square();
        let t2 = self.z.square();
        let t3 = &self.x * &self.y;
        let t3 = t3.double();
        let z3 = &self.x * &self.z;
        let z3 = &z3 + &z3;
        let x3 = curve.a() * &z3;
        let y3 = curve.b3() * &t2;
        let y3 = &x3 + &y3;
        let x3 = &t1 - &y3;
        let y3 = &t1 + &y3;
        let y3 = &x3 * &y3;
        let x3 = &t3 * &x3;
        let z3 = curve.b3() * &z3;
        let t2 = curve.a() * &t2;
        let t3 = &t0 - &t2;
        let t3 = curve.a() * &t3;
        let t3 = &t3 + &z3;
        let z3 = &t0 + &t0;
        let t0 = &z3 + &t0;
        let t0 = &t0 + &t2;
        let t0 = &t0 * &t3;
        let y3 = &y3 + &t0;
        let t2 = &self.y * &self.z;
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

    #[inline]
    pub fn double_a0<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveA0>(
        &self,
        curve: C,
    ) -> Self {
        // Algorithm 5 from (1) - doubling formula for a=0
        //
        // ```magma
        // DBL := function (X ,Y ,Z , b3 )
        //     t0 := Y ^2; Z3 := t0 + t0 ; Z3 := Z3 + Z3 ;
        //     Z3 := Z3 + Z3 ; t1 := Y * Z ; t2 := Z ^2;
        //     t2 := b3 * t2 ; X3 := t2 * Z3 ; Y3 := t0 + t2 ;
        //     Z3 := t1 * Z3 ; t1 := t2 + t2 ; t2 := t1 + t2 ;
        //     t0 := t0 - t2 ; Y3 := t0 * Y3 ; Y3 := X3 + Y3 ;
        //     t1 := X * Y ; X3 := t0 * t1 ; X3 := X3 + X3 ;
        //     return X3 , Y3 , Z3 ;
        // end function ;
        // ```

        let t0 = self.y.square();
        let z3 = &t0 + &t0;
        let z3 = z3.double();
        let z3 = z3.double();
        let t1 = &self.y * &self.z;
        let t2 = self.z.square();
        let t2 = curve.b3() * &t2;
        let x3 = &t2 * &z3;
        let y3 = &t0 + &t2;
        let z3 = &t1 * &z3;
        let t1 = &t2 + &t2;
        let t2 = &t1 + &t2;
        let t0 = &t0 - &t2;
        let y3 = &t0 * &y3;
        let y3 = &x3 + &y3;
        let t1 = &self.x * &self.y;
        let x3 = &t0 * &t1;
        let x3 = x3.double();

        Point {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    pub fn to_affine(&self) -> Option<affine::Point<FE>> {
        if self.z == FE::one() {
            return Some(affine::Point {
                x: self.x.clone(),
                y: self.y.clone(),
            });
        }
        if self.z.is_zero() {
            None
        } else {
            let inv = self.z.inverse();
            Some(affine::Point {
                x: &self.x * &inv,
                y: &self.y * &inv,
            })
        }
    }

    /// scalar multiplication : `n * self` with double-and-add algorithm with increasing index
    #[inline]
    fn scalar_mul_daa_limbs8<C: WeierstrassCurve<FieldElement = FE>>(
        &self,
        n: &[u8],
        curve: C,
    ) -> Self {
        let mut a: Point<FE> = self.clone();
        let mut q: Point<FE> = Point::infinity();

        for digit in n.iter().rev() {
            for i in 0..8 {
                if digit & (1 << i) != 0 {
                    q = q.add_or_double(&a, curve);
                }
                a = a.double(curve)
            }
        }
        q
    }

    #[inline]
    fn scalar_mul_daa_limbs8_a0<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveA0>(
        &self,
        n: &[u8],
        curve: C,
    ) -> Self {
        let mut a: Point<FE> = self.clone();
        let mut q: Point<FE> = Point::infinity();

        for digit in n.iter().rev() {
            for i in 0..8 {
                if digit & (1 << i) != 0 {
                    q = q.add_or_double_a0(&a, curve);
                }
                a = a.double_a0(curve)
            }
        }
        q
    }

    pub fn scale<C: WeierstrassCurve<FieldElement = FE>>(&self, n: &[u8], curve: C) -> Self {
        self.scalar_mul_daa_limbs8(n, curve)
    }

    pub fn scale_a0<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveA0>(
        &self,
        n: &[u8],
        curve: C,
    ) -> Self {
        self.scalar_mul_daa_limbs8_a0(n, curve)
    }

    #[inline]
    pub fn add_or_double<'b, C: WeierstrassCurve<FieldElement = FE>>(
        &self,
        other: &'b Point<FE>,
        curve: C,
    ) -> Point<FE> {
        if self.ct_eq(other).is_true() {
            self.double(curve)
        } else {
            self.add_different(&other, curve)
        }
    }

    #[inline]
    pub fn add_or_double_a0<'b, C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveA0>(
        &self,
        other: &'b Point<FE>,
        curve: C,
    ) -> Point<FE> {
        if self.ct_eq(other).is_true() {
            self.double_a0(curve)
        } else {
            self.add_different_a0(&other, curve)
        }
    }
}

impl<FE> std::ops::Neg for Point<FE>
where
    FE: Neg<Output = FE>,
{
    type Output = Point<FE>;

    fn neg(self) -> Self::Output {
        Point {
            x: self.x,
            y: -self.y,
            z: self.z,
        }
    }
}

impl<'a, FE> std::ops::Neg for &'a Point<FE>
where
    FE: Clone + Neg<Output = FE>,
    &'a FE: Neg<Output = FE>,
{
    type Output = Point<FE>;

    fn neg(self) -> Self::Output {
        Point {
            x: self.x.clone(),
            y: -&self.y,
            z: self.z.clone(),
        }
    }
}
