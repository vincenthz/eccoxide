//! Jacobian Elliptic Curve Point defined over Field element as (X,Y,Z)
//!
//! A jacobian point `(X:Y:Z)` with `Z != 0` is the affine point
//! `(X/Z², Y/Z³)`, and `Z == 0` is the point at infinity.
//!
//! Compared with the complete homogeneous projective formulas of
//! [`super::projective`]:
//!
//! * the doubling is markedly cheaper; `3M + 5S` for `a = -3`
//!   and `2M + 5S` for `a = 0`, against `8M + 3S + 2·m_b`
//! * an addend already in affine form (`Z2 = 1`) gets a dedicated mixed
//!   addition (`7M + 4S`) that has no equivalent with the complete formulas.
//!
//! However the formulas are incomplete and give wrong result when:
//!
//! * The two inputs are the same point
//! * Either input are at infinity
//!
//! As a results, the APIs available:
//!
//! * [`Point::add_different`] / [`Point::add_different_mixed`]: the raw
//!   formulas, for callers which know the exceptional cases cannot happen
//! * [`Point::add_or_double`]: constant-time wrappers handling every
//!   case, at the cost of also computing the unnecessary operation
//!
//! The formulas are the ones of the
//! [Explicit-Formulas Database](https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html):
//! `add-2007-bl`, `madd-2007-bl`, `dbl-2007-bl` (arbitrary `a`), `dbl-2009-l`
//! (`a = 0`) and `dbl-2001-b` (`a = -3`).
//!

use super::affine;
use super::field::Field;
use super::projective;
use super::weierstrass::{WeierstrassCurve, WeierstrassCurveA0, WeierstrassCurveAM3};
use super::wnaf::wnaf;
use crate::mp::ct::{Choice, CtEqual, CtOption, CtSelect};
use alloc::vec::Vec;
use core::convert::TryFrom;
use core::ops::{Add, Mul, Neg, Sub};

pub use super::projective::AffineAtInfinity;

/// Jacobian point with field element FE
///
/// Affine point associated with (X,Y,Z) : (X/Z², Y/Z³)
///
/// Note that 2 points are equal if they are in the same equivalence class,
/// which is determined with 6 FieldElement multiplications.
///
/// Example: (1,1,1) and (4,8,2) are equal
#[derive(Clone, Debug)]
pub struct Point<FE> {
    pub x: FE,
    pub y: FE,
    pub z: FE,
}

impl<FE: Field> PartialEq for Point<FE>
where
    for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
{
    fn eq(&self, other: &Point<FE>) -> bool {
        self.is_equivalent(other).is_true()
    }
}

impl<FE: Field> Eq for Point<FE> where for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE> {}

impl<FE: Field> CtEqual for Point<FE>
where
    for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
{
    fn ct_eq(&self, other: &Point<FE>) -> Choice {
        self.is_equivalent(other)
    }
}

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
    FE: Field,
    for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
{
    /// Check if two points are in the same equivalence class
    ///
    /// `(X1:Y1:Z1)` and `(X2:Y2:Z2)` are the same affine point when
    /// `X1*Z2² = X2*Z1²` and `Y1*Z2³ = Y2*Z1³`. Those two conditions alone
    /// also hold between the degenerate `(0:0:0)` and any point of
    /// x-coordinate zero, so the infinity flags are compared as well.
    fn is_equivalent(&self, other: &Point<FE>) -> Choice {
        let z1z1 = self.z.square();
        let z2z2 = other.z.square();
        let nx1 = &self.x * &z2z2;
        let nx2 = &other.x * &z1z1;
        let ny1 = &self.y * &(&z2z2 * &other.z);
        let ny2 = &other.y * &(&z1z1 * &self.z);
        nx1.ct_eq(&nx2) & ny1.ct_eq(&ny2) & self.is_infinity().ct_eq(&other.is_infinity())
    }
}

impl<FE> Point<FE>
where
    FE: Field,
{
    /// The point at infinity
    pub const INFINITY: Self = Point {
        x: FE::ONE,
        y: FE::ONE,
        z: FE::ZERO,
    };

    /// Check if a point is at infinity
    pub fn is_infinity(&self) -> Choice {
        self.z.ct_eq(&FE::ZERO)
    }

    pub fn from_affine(p: &affine::Point<FE>) -> Self {
        Point {
            x: p.x.clone(),
            y: p.y.clone(),
            z: FE::ONE,
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
            let zinv2 = zinv.square();
            let zinv3 = &zinv2 * &zinv;

            self.x = &self.x * &zinv2;
            self.y = &self.y * &zinv3;
            self.z = FE::ONE
        }
    }

    /// Convert a homogeneous projective point to its jacobian equivalent
    ///
    /// `(X:Y:Z)` projective is the affine point `(X/Z, Y/Z)`, which as a
    /// jacobian point of same `Z` is `(X*Z : Y*Z² : Z)`.
    pub fn from_projective(p: &projective::Point<FE>) -> Self {
        let q = Point {
            x: &p.x * &p.z,
            y: &p.y * &p.z.square(),
            z: p.z.clone(),
        };
        // the projective point at infinity has X = 0, so the conversion of it
        // is the degenerate (0:0:0) rather than INFINITY
        Self::ct_select(p.z.ct_eq(&FE::ZERO), &Self::INFINITY, &q)
    }

    /// Convert this point to its homogeneous projective equivalent
    ///
    /// `(X:Y:Z)` jacobian is the affine point `(X/Z², Y/Z³)`, which as a
    /// projective point of `Z³` is `(X·Z : Y : Z³)`. The point at infinity
    /// maps to the projective `(0:1:0)`.
    pub fn to_projective(&self) -> projective::Point<FE> {
        let zz = self.z.square();
        projective::Point {
            x: &self.x * &self.z,
            y: self.y.clone(),
            z: &zz * &self.z,
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
    /// Add two *distinct* points, neither of them at infinity: `11M + 5S`
    ///
    /// The result is not the sum when `self` and `other` are the same point
    /// (the formula degenerates to `(0:0:0)`) or when either of them is at
    /// infinity. Use [`Self::add_or_double`] when either can happen.
    #[inline]
    pub fn add_different(&self, other: &Point<FE>) -> Point<FE> {
        self.add_inner(other).0
    }

    /// Add a *distinct* point given in affine form (`other.z == 1`), neither
    /// of them at infinity: `7M + 4S`
    ///
    /// Same exceptional cases as [`Self::add_different`]; see
    /// [`Self::add_or_double_mixed`] for the complete variant.
    #[inline]
    pub fn add_different_mixed(&self, other: &Point<FE>) -> Point<FE> {
        self.add_inner_mixed(other).0
    }

    /// Addition formula `add-2007-bl` from (1), returning the sum along with
    /// whether the two inputs turned out to be the same point.
    ///
    /// ```text
    /// ADD := function (X1, Y1, Z1, X2, Y2, Z2)
    ///     Z1Z1 := Z1^2;
    ///     Z2Z2 := Z2^2;
    ///     U1 := X1 * Z2Z2;
    ///     U2 := X2 * Z1Z1;
    ///     S1 := Y1 * Z2 * Z2Z2;
    ///     S2 := Y2 * Z1 * Z1Z1;
    ///     H := U2 - U1;
    ///     I := (2*H)^2;
    ///     J := H * I;
    ///     r := 2*(S2 - S1);
    ///     V := U1 * I;
    ///     X3 := r^2 - J - 2*V;
    ///     Y3 := r*(V - X3) - 2*S1*J;
    ///     Z3 := ((Z1 + Z2)^2 - Z1Z1 - Z2Z2) * H;
    ///     return X3, Y3, Z3;
    /// end function;
    /// ```
    ///
    /// `H = 0` means the two points share their affine x-coordinate, and
    /// `r = 0` on top of it that they share their affine y-coordinate too:
    /// together they are exactly the doubling case, which the formula cannot
    /// compute (every output coordinate is zero).
    fn add_inner(&self, other: &Point<FE>) -> (Point<FE>, Choice) {
        let z1z1 = self.z.square();
        let z2z2 = other.z.square();
        let u1 = &self.x * &z2z2;
        let u2 = &other.x * &z1z1;
        let s1 = &self.y * &other.z;
        let s1 = &s1 * &z2z2;
        let s2 = &other.y * &self.z;
        let s2 = &s2 * &z1z1;
        let h = &u2 - &u1;
        let i = h.double().square();
        let j = &h * &i;
        let r = (&s2 - &s1).double();
        let v = &u1 * &i;

        let x3 = &r.square() - &j;
        let x3 = &x3 - &v.double();
        let y3 = &r * &(&v - &x3);
        let y3 = &y3 - &(&s1 * &j).double();
        let z3 = &(&self.z + &other.z).square() - &z1z1;
        let z3 = &z3 - &z2z2;
        let z3 = &z3 * &h;

        let same = h.ct_eq(&FE::ZERO) & r.ct_eq(&FE::ZERO);
        (
            Point {
                x: x3,
                y: y3,
                z: z3,
            },
            same,
        )
    }

    /// Mixed addition formula `madd-2007-bl` from (1) — `other` is read as the
    /// affine point `(X2, Y2)`, i.e. its `z` is assumed to be one — returning
    /// the sum along with whether the two inputs are the same point.
    ///
    /// ```text
    /// MADD := function (X1, Y1, Z1, X2, Y2)
    ///     Z1Z1 := Z1^2;
    ///     U2 := X2 * Z1Z1;
    ///     S2 := Y2 * Z1 * Z1Z1;
    ///     H := U2 - X1;
    ///     HH := H^2;
    ///     I := 4*HH;
    ///     J := H * I;
    ///     r := 2*(S2 - Y1);
    ///     V := X1 * I;
    ///     X3 := r^2 - J - 2*V;
    ///     Y3 := r*(V - X3) - 2*Y1*J;
    ///     Z3 := (Z1 + H)^2 - Z1Z1 - HH;
    ///     return X3, Y3, Z3;
    /// end function;
    /// ```
    fn add_inner_mixed(&self, other: &Point<FE>) -> (Point<FE>, Choice) {
        let z1z1 = self.z.square();
        let u2 = &other.x * &z1z1;
        let s2 = &other.y * &self.z;
        let s2 = &s2 * &z1z1;
        let h = &u2 - &self.x;
        let hh = h.square();
        let i = hh.double().double();
        let j = &h * &i;
        let r = (&s2 - &self.y).double();
        let v = &self.x * &i;

        let x3 = &r.square() - &j;
        let x3 = &x3 - &v.double();
        let y3 = &r * &(&v - &x3);
        let y3 = &y3 - &(&self.y * &j).double();
        let z3 = &(&self.z + &h).square() - &z1z1;
        let z3 = &z3 - &hh;

        let same = h.ct_eq(&FE::ZERO) & r.ct_eq(&FE::ZERO);
        (
            Point {
                x: x3,
                y: y3,
                z: z3,
            },
            same,
        )
    }

    /// Turn an incomplete addition into a complete one, in constant time.
    #[inline]
    fn complete<D>(&self, other: &Point<FE>, sum: (Point<FE>, Choice), dbl: D) -> Point<FE>
    where
        D: FnOnce(&Point<FE>) -> Point<FE>,
    {
        let (sum, same) = sum;
        let p = Self::ct_select(same, &dbl(self), &sum);
        let p = Self::ct_select(other.is_infinity(), self, &p);
        Self::ct_select(self.is_infinity(), other, &p)
    }

    /// Variable-time counterpart of [`Self::complete`]: the exceptional cases
    /// are branched on rather than selected, so neither the addition nor the
    /// doubling is evaluated unless it is the answer.
    ///
    /// Only for public point values
    #[inline]
    fn add_or_double_vartime<D>(&self, other: &Point<FE>, dbl: D) -> Point<FE>
    where
        D: FnOnce(&Point<FE>) -> Point<FE>,
    {
        if self.is_infinity().is_true() {
            return other.clone();
        }
        if other.is_infinity().is_true() {
            return self.clone();
        }
        let (sum, same) = self.add_inner(other);
        if same.is_true() {
            dbl(self)
        } else {
            sum
        }
    }

    /// Add two points, correctly handling every case (`self == other`,
    /// `self == -other` and the point at infinity), in constant time.
    ///
    /// Unlike the complete formulas of [`super::projective`] this is not a
    /// single exception-free formula: it evaluates the incomplete addition
    /// *and* the doubling and selects between them, so it costs more than
    /// either. Prefer [`Self::add_different`] wherever the operands are known
    /// to be distinct and finite.
    #[inline]
    pub fn add_or_double<C: WeierstrassCurve<FieldElement = FE>>(
        &self,
        other: &Point<FE>,
    ) -> Point<FE> {
        self.complete(other, self.add_inner(other), |p| p.double::<C>())
    }

    /// [`Self::add_or_double`] for a=0 curves
    #[inline]
    pub fn add_or_double_a0<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveA0>(
        &self,
        other: &Point<FE>,
    ) -> Point<FE> {
        self.complete(other, self.add_inner(other), |p| p.double_a0::<C>())
    }

    /// [`Self::add_or_double`] for a=-3 curves
    #[inline]
    pub fn add_or_double_am3<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveAM3>(
        &self,
        other: &Point<FE>,
    ) -> Point<FE> {
        self.complete(other, self.add_inner(other), |p| p.double_am3::<C>())
    }

    /// [`Self::add_or_double`] with `other` in affine form (`other.z == 1`),
    /// or at infinity.
    #[inline]
    pub fn add_or_double_mixed<C: WeierstrassCurve<FieldElement = FE>>(
        &self,
        other: &Point<FE>,
    ) -> Point<FE> {
        self.complete(other, self.add_inner_mixed(other), |p| p.double::<C>())
    }

    /// [`Self::add_or_double_mixed`] for a=0 curves
    #[inline]
    pub fn add_or_double_mixed_a0<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveA0>(
        &self,
        other: &Point<FE>,
    ) -> Point<FE> {
        self.complete(other, self.add_inner_mixed(other), |p| p.double_a0::<C>())
    }

    /// [`Self::add_or_double_mixed`] for a=-3 curves
    #[inline]
    pub fn add_or_double_mixed_am3<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveAM3>(
        &self,
        other: &Point<FE>,
    ) -> Point<FE> {
        self.complete(other, self.add_inner_mixed(other), |p| p.double_am3::<C>())
    }

    #[inline]
    pub fn double<C: WeierstrassCurve<FieldElement = FE>>(&self) -> Self {
        // Algorithm `dbl-2007-bl` from (1) - doubling for arbitrary a : 1M + 8S + 1*a
        //
        // ```text
        // DBL := function (X, Y, Z, a)
        //     XX := X^2;
        //     YY := Y^2;
        //     YYYY := YY^2;
        //     ZZ := Z^2;
        //     S := 2*((X + YY)^2 - XX - YYYY);
        //     M := 3*XX + a*ZZ^2;
        //     T := M^2 - 2*S;
        //     X3 := T;
        //     Y3 := M*(S - T) - 8*YYYY;
        //     Z3 := (Y + Z)^2 - YY - ZZ;
        //     return X3, Y3, Z3;
        // end function;
        // ```
        let xx = self.x.square();
        let yy = self.y.square();
        let yyyy = yy.square();
        let zz = self.z.square();
        let s = &(&self.x + &yy).square() - &xx;
        let s = (&s - &yyyy).double();
        let m = &xx.double() + &xx;
        let m = &m + &(C::A * &zz.square());
        let t = &m.square() - &s.double();
        let y3 = &(&m * &(&s - &t)) - &yyyy.double().double().double();
        let z3 = &(&self.y + &self.z).square() - &yy;
        let z3 = &z3 - &zz;

        Point { x: t, y: y3, z: z3 }
    }

    #[inline]
    pub fn double_a0<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveA0>(&self) -> Self {
        // Algorithm `dbl-2009-l` from (1) - doubling for a=0 : 2M + 5S
        //
        // ```text
        // DBL := function (X, Y, Z)
        //     A := X^2;
        //     B := Y^2;
        //     C := B^2;
        //     D := 2*((X + B)^2 - A - C);
        //     E := 3*A ; F := E^2;
        //     X3 := F - 2*D;
        //     Y3 := E*(D - X3) - 8*C;
        //     Z3 := 2*Y*Z;
        //     return X3, Y3, Z3;
        // end function;
        // ```
        let a = self.x.square();
        let b = self.y.square();
        let c = b.square();
        let d = &(&self.x + &b).square() - &a;
        let d = (&d - &c).double();
        let e = &a.double() + &a;
        let f = e.square();
        let x3 = &f - &d.double();
        let y3 = &(&e * &(&d - &x3)) - &c.double().double().double();
        let z3 = (&self.y * &self.z).double();

        Point {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    #[inline]
    pub fn double_am3<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveAM3>(&self) -> Self {
        // Algorithm `dbl-2001-b` from (1) - doubling for a=-3 : 3M + 5S
        //
        // ```text
        // DBL := function (X, Y, Z)
        //     delta := Z^2;
        //     gamma := Y^2;
        //     beta := X * gamma;
        //     alpha := 3*(X - delta)*(X + delta);
        //     X3 := alpha^2 - 8*beta;
        //     Z3 := (Y + Z)^2 - gamma - delta;
        //     Y3 := alpha*(4*beta - X3) - 8*gamma^2;
        //     return X3, Y3, Z3;
        // end function;
        // ```
        let delta = self.z.square();
        let gamma = self.y.square();
        let beta = &self.x * &gamma;
        let alpha = &(&self.x - &delta) * &(&self.x + &delta);
        let alpha = &alpha.double() + &alpha;
        let beta4 = beta.double().double();
        let x3 = &alpha.square() - &beta4.double();
        let z3 = &(&self.y + &self.z).square() - &gamma;
        let z3 = &z3 - &delta;
        let y3 = &(&alpha * &(&beta4 - &x3)) - &gamma.square().double().double().double();

        Point {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    /// Presence choice `z != 0` together with `1/z`, `1/z²` and `1/z³`.
    ///
    /// For the point at infinity (`z == 0`) the inversion input is
    /// substituted with 1 so that it stays defined; the returned inverses are
    /// then placeholders to be discarded through the false choice. Shared by
    /// the constant-time affine conversions.
    #[inline]
    fn z_inverse_ct(&self) -> (Choice, FE, FE) {
        let z_nonzero = self.z.ct_ne(&FE::ZERO);
        let z = FE::ct_select(z_nonzero, &self.z, &FE::ONE);
        let zinv = z.inverse();
        let zinv2 = zinv.square();
        let zinv3 = &zinv2 * &zinv;
        (z_nonzero, zinv2, zinv3)
    }

    /// Try to return the affine Point for this jacobian point
    ///
    /// If Self represent the point at infinity then None is return.
    ///
    /// For a constant-time variant of this, use [`Self::to_affine_ct`]
    pub fn to_affine(&self) -> Option<affine::Point<FE>> {
        self.to_affine_ct().into_option()
    }

    /// Constant-time variant of [`Self::to_affine`]: no branch on the point
    /// value.
    ///
    /// The presence choice is `z != 0`; for the point at infinity the carried
    /// coordinates are placeholders to be discarded (see [`Self::z_inverse_ct`]).
    pub fn to_affine_ct(&self) -> CtOption<affine::Point<FE>> {
        let (z_nonzero, zinv2, zinv3) = self.z_inverse_ct();
        let p = affine::Point {
            x: &self.x * &zinv2,
            y: &self.y * &zinv3,
        };
        CtOption::from((z_nonzero, p))
    }

    /// Constant-time affine x-coordinate: like [`Self::to_affine_ct`] but
    /// computes only `x/z²`, skipping the y-coordinate when it is not needed.
    ///
    /// The presence choice is `z != 0`; for the point at infinity the carried
    /// x is a placeholder to be discarded (see [`Self::z_inverse_ct`]).
    pub fn to_affine_x_ct(&self) -> CtOption<FE> {
        let (z_nonzero, zinv2, _) = self.z_inverse_ct();
        CtOption::from((z_nonzero, &self.x * &zinv2))
    }

    /// Additive inverse of the point: `-(X:Y:Z) = (X:-Y:Z)`.
    #[inline]
    fn negate(&self) -> Self {
        Point {
            x: self.x.clone(),
            y: -self.y.clone(),
            z: self.z.clone(),
        }
    }

    /// Constant-time lookup of `table[index]`. The whole table is scanned so the
    /// memory access pattern does not depend on the (secret) `index`.
    fn select_from_table(table: &[Point<FE>], index: u8) -> Point<FE> {
        let mut acc = Self::INFINITY;
        for (j, t) in table.iter().enumerate() {
            let take = (j as u64).ct_eq(&(index as u64));
            acc = Point::ct_select(take, t, &acc);
        }
        acc
    }

    /// Width of the signed window used by the variable-time wNAF scalar
    /// multiplication. `W = 5` keeps the odd-multiple table small (`2^(W-2)`
    /// points) while giving an average non-zero digit density of `1/(W+1)`,
    /// which is close to optimal for ~192–521 bit scalars.
    const WNAF_W: u32 = 5;

    /// Variable-time scalar multiplication `n * self` using a width-`W` wNAF.
    ///
    /// Compared with a plain double-and-add this keeps the same number of
    /// doublings but cuts the additions from ~`nbits/2` down to ~`nbits/(W+1)`
    /// by using a signed sparse recoding of the scalar. The running time
    /// depends on the scalar, so this must only be used when `n` is public.
    ///
    /// The doubling and the (variable-time complete) addition are taken as
    /// parameters so the a=0 / a=-3 / arbitrary-a variants share this body.
    #[inline]
    fn scalar_mul_wnaf<D, A>(&self, n: &[u8], dbl: D, add: A) -> Self
    where
        D: Fn(&Self) -> Self,
        A: Fn(&Self, &Self) -> Self,
    {
        let naf = wnaf(n, Self::WNAF_W);
        // table[i] = (2i+1) * self, i.e. the odd multiples 1·P, 3·P, … of self
        let dbl2 = dbl(self);
        let tlen = 1usize << (Self::WNAF_W - 2);
        let mut table: Vec<Point<FE>> = Vec::with_capacity(tlen);
        table.push(self.clone());
        for i in 1..tlen {
            table.push(add(&table[i - 1], &dbl2));
        }

        let mut q = Self::INFINITY;
        for &d in naf.iter().rev() {
            q = dbl(&q);
            if d > 0 {
                q = add(&q, &table[(d as usize) >> 1]);
            } else if d < 0 {
                q = add(&q, &table[((-(d as i32)) as usize) >> 1].negate());
            }
        }
        q
    }

    /// Constant-time scalar multiplication using a fixed 4-bit window.
    ///
    /// The window value is read from a precomputed table with a constant-time
    /// scan, and the addition is the complete wrapper so that a zero window,
    /// a zero accumulator or an accumulator that happens to equal the window
    /// point all give the right answer without branching.
    ///
    /// The doubling and the addition are taken as parameters so the a=0 /
    /// a=-3 / arbitrary-a variants share this body.
    fn scalar_mul_fixed_window<D, A>(&self, n: &[u8], dbl: D, add: A) -> Self
    where
        D: Fn(&Self) -> Self,
        A: Fn(&Self, &Self) -> Self,
    {
        // table[d] = d * self, for d in 0..16
        let mut table: [Point<FE>; 16] = core::array::from_fn(|_| Self::INFINITY);
        table[1] = self.clone();
        table[2] = dbl(self);
        for d in 3..16 {
            let next = add(&table[d - 1], self);
            table[d] = next;
        }

        let mut q = Self::INFINITY;
        for byte in n.iter() {
            for &index in &[byte >> 4, byte & 0x0f] {
                q = dbl(&dbl(&dbl(&dbl(&q))));
                let selected = Self::select_from_table(&table, index);
                q = add(&q, &selected);
            }
        }
        q
    }

    /// Variable-time scalar multiplication for arbitrary a curves.
    /// See [`Self::scalar_mul_wnaf`].
    pub fn scale<C: WeierstrassCurve<FieldElement = FE>>(&self, n: &[u8]) -> Self {
        self.scalar_mul_wnaf(
            n,
            |p| p.double::<C>(),
            |p, q| p.add_or_double_vartime(q, |p| p.double::<C>()),
        )
    }

    /// Variable-time scalar multiplication for a=0 curves.
    pub fn scale_a0<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveA0>(
        &self,
        n: &[u8],
    ) -> Self {
        self.scalar_mul_wnaf(
            n,
            |p| p.double_a0::<C>(),
            |p, q| p.add_or_double_vartime(q, |p| p.double_a0::<C>()),
        )
    }

    /// Variable-time scalar multiplication for a=-3 curves.
    pub fn scale_am3<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveAM3>(
        &self,
        n: &[u8],
    ) -> Self {
        self.scalar_mul_wnaf(
            n,
            |p| p.double_am3::<C>(),
            |p, q| p.add_or_double_vartime(q, |p| p.double_am3::<C>()),
        )
    }

    /// Constant-time scalar multiplication (default). See
    /// [`Self::scalar_mul_fixed_window`].
    pub fn scale_ct<C: WeierstrassCurve<FieldElement = FE>>(&self, n: &[u8]) -> Self {
        self.scalar_mul_fixed_window(n, |p| p.double::<C>(), |p, q| p.add_or_double::<C>(q))
    }

    /// Constant-time scalar multiplication for a=0 curves (default).
    pub fn scale_a0_ct<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveA0>(
        &self,
        n: &[u8],
    ) -> Self {
        self.scalar_mul_fixed_window(n, |p| p.double_a0::<C>(), |p, q| p.add_or_double_a0::<C>(q))
    }

    /// Constant-time scalar multiplication for a=-3 curves (default).
    pub fn scale_am3_ct<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveAM3>(
        &self,
        n: &[u8],
    ) -> Self {
        self.scalar_mul_fixed_window(
            n,
            |p| p.double_am3::<C>(),
            |p, q| p.add_or_double_am3::<C>(q),
        )
    }
}

impl<FE: Field> CtSelect for Point<FE> {
    /// Constant-time select between two points: returns `a` if `cond` is true,
    /// otherwise `b`, without branching on `cond`.
    fn ct_select(cond: Choice, a: &Point<FE>, b: &Point<FE>) -> Point<FE> {
        Point {
            x: FE::ct_select(cond, &a.x, &b.x),
            y: FE::ct_select(cond, &a.y, &b.y),
            z: FE::ct_select(cond, &a.z, &b.z),
        }
    }
    fn ct_assign(&mut self, cond: Choice, other: &Point<FE>) {
        self.x.ct_assign(cond, &other.x);
        self.y.ct_assign(cond, &other.y);
        self.z.ct_assign(cond, &other.z);
    }
}

impl<FE> core::ops::Neg for Point<FE>
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

impl<'a, FE> core::ops::Neg for &'a Point<FE>
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Exercise the representation and the arbitrary-`a` formulas of a curve,
    /// checking everything against the (complete, independently derived)
    /// projective implementation of the same curve.
    fn check_curve<FE, C>(g: &affine::Point<FE>)
    where
        FE: Field,
        C: WeierstrassCurve<FieldElement = FE>,
        for<'a> &'a FE: Add<FE, Output = FE>,
        for<'a> &'a FE: Mul<FE, Output = FE>,
        for<'a> &'a FE: Sub<FE, Output = FE>,
        for<'a, 'b> &'a FE: Add<&'b FE, Output = FE>,
        for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
        for<'a, 'b> &'a FE: Sub<&'b FE, Output = FE>,
    {
        let inf = Point::<FE>::INFINITY;
        let g1 = Point::from_affine(g);
        let p1 = projective::Point::from_affine(g);

        // affine conversion round trip, and the point at infinity having none
        assert_eq!(g1.to_affine().as_ref(), Some(g));
        assert!(g1.is_infinity().is_false());
        assert!(inf.is_infinity().is_true());
        assert!(inf.to_affine().is_none());
        assert!(inf.to_affine_ct().is_present().is_false());
        assert_eq!(g1.to_affine_x_ct().into_option().as_ref(), Some(&g.x));

        // doubling and addition against the projective implementation
        let g2 = g1.double::<C>();
        let q2 = p1.double::<C>();
        assert_eq!(g2.to_affine(), q2.to_affine());

        let g3 = g2.add_different(&g1);
        let q3 = q2.add_different::<C>(&p1);
        assert_eq!(g3.to_affine(), q3.to_affine());

        // the mixed addition agrees with the general one, `g1` having z = 1
        assert_eq!(g2.add_different_mixed(&g1).to_affine(), g3.to_affine());

        // 4G through two different routes
        let g4 = g2.double::<C>();
        assert_eq!(g4.to_affine(), g3.add_different(&g1).to_affine());

        // every exceptional case of the incomplete formulas, through the
        // complete wrappers, for both the general and the mixed addition
        for &mixed in &[false, true] {
            let add = |a: &Point<FE>, b: &Point<FE>| {
                if mixed {
                    a.add_or_double_mixed::<C>(b)
                } else {
                    a.add_or_double::<C>(b)
                }
            };
            // P + P is the doubling
            assert_eq!(add(&g1, &g1).to_affine(), g2.to_affine());
            // P + (-P) is the point at infinity
            assert!(add(&g1, &(-g1.clone())).is_infinity().is_true());
            // the point at infinity is the neutral element on both sides
            assert_eq!(add(&g1, &inf).to_affine(), g1.to_affine());
            assert_eq!(add(&inf, &g1).to_affine(), g1.to_affine());
            assert!(add(&inf, &inf).is_infinity().is_true());
            // and the ordinary case still is the ordinary case
            assert_eq!(add(&g2, &g1).to_affine(), g3.to_affine());
        }

        // a point and a rescaled representative of it are equal, and
        // normalizing one gives the z = 1 representative
        let l = FE::from(7u64);
        let l2 = l.square();
        let l3 = &l2 * &l;
        let rescaled = Point {
            x: &g3.x * &l2,
            y: &g3.y * &l3,
            z: &g3.z * &l,
        };
        assert_eq!(rescaled, g3);
        assert!(rescaled.ct_eq(&g3).is_true());
        assert_ne!(rescaled, g4);
        let mut normalized = rescaled.clone();
        normalized.normalize();
        assert_eq!(normalized.z, FE::ONE);
        assert_eq!(normalized, g3);
        assert_eq!(normalized.to_affine(), g3.to_affine());

        // the degenerate (0:0:0) the incomplete addition returns for a
        // doubling is at infinity, and must not compare equal to a finite
        // point of x-coordinate zero — which is what the coordinate products
        // alone would say
        let zero = Point {
            x: FE::ZERO,
            y: FE::ZERO,
            z: FE::ZERO,
        };
        let x_zero = Point {
            x: FE::ZERO,
            y: FE::ONE,
            z: FE::ONE,
        };
        assert!(zero.is_infinity().is_true());
        assert_eq!(zero, inf);
        assert_ne!(zero, x_zero);
        assert_ne!(x_zero, zero);
        assert_ne!(inf, g1);
        assert_ne!(g1, inf);

        // interoperability with the projective representation
        assert_eq!(Point::from_projective(&q3).to_affine(), g3.to_affine());
        assert_eq!(g3.to_projective().to_affine(), q3.to_affine());
        assert!(Point::<FE>::from_projective(&projective::Point::INFINITY)
            .is_infinity()
            .is_true());
        assert!(inf.to_projective().is_infinity().is_true());

        // scalar multiplication, constant-time and variable-time, against the
        // projective implementation and against repeated addition
        assert!(g1.scale_ct::<C>(&[0]).is_infinity().is_true());
        assert!(g1.scale::<C>(&[0]).is_infinity().is_true());
        assert_eq!(g1.scale_ct::<C>(&[1]).to_affine(), g1.to_affine());
        assert_eq!(g1.scale::<C>(&[1]).to_affine(), g1.to_affine());
        assert_eq!(g1.scale_ct::<C>(&[4]).to_affine(), g4.to_affine());
        assert_eq!(g1.scale::<C>(&[4]).to_affine(), g4.to_affine());

        let n: [u8; 8] = [0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xf0];
        let expected = p1.scale::<C>(&n).to_affine();
        assert_eq!(p1.scale_ct::<C>(&n).to_affine(), expected);
        assert_eq!(g1.scale::<C>(&n).to_affine(), expected);
        assert_eq!(g1.scale_ct::<C>(&n).to_affine(), expected);

        // scalar multiplication of the point at infinity stays there
        assert!(inf.scale_ct::<C>(&n).is_infinity().is_true());
        assert!(inf.scale::<C>(&n).is_infinity().is_true());
    }

    /// The a=0 specialised formulas must agree with the arbitrary-`a` ones
    fn check_a0<FE, C>(g: &affine::Point<FE>)
    where
        FE: Field,
        C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveA0,
        for<'a> &'a FE: Add<FE, Output = FE>,
        for<'a> &'a FE: Mul<FE, Output = FE>,
        for<'a> &'a FE: Sub<FE, Output = FE>,
        for<'a, 'b> &'a FE: Add<&'b FE, Output = FE>,
        for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
        for<'a, 'b> &'a FE: Sub<&'b FE, Output = FE>,
    {
        let inf = Point::<FE>::INFINITY;
        let g1 = Point::from_affine(g);
        let g2 = g1.double::<C>();

        assert_eq!(g1.double_a0::<C>().to_affine(), g2.to_affine());
        assert_eq!(
            g2.double_a0::<C>().to_affine(),
            g2.double::<C>().to_affine()
        );
        assert!(inf.double_a0::<C>().is_infinity().is_true());

        assert_eq!(
            g2.add_or_double_a0::<C>(&g1).to_affine(),
            g2.add_or_double::<C>(&g1).to_affine()
        );
        assert_eq!(g1.add_or_double_a0::<C>(&g1).to_affine(), g2.to_affine());
        assert_eq!(
            g2.add_or_double_mixed_a0::<C>(&g1).to_affine(),
            g2.add_or_double::<C>(&g1).to_affine()
        );
        assert!(g1
            .add_or_double_a0::<C>(&(-g1.clone()))
            .is_infinity()
            .is_true());

        let n: [u8; 8] = [0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xf0];
        let expected = g1.scale_ct::<C>(&n).to_affine();
        assert_eq!(g1.scale_a0_ct::<C>(&n).to_affine(), expected);
        assert_eq!(g1.scale_a0::<C>(&n).to_affine(), expected);
    }

    /// The a=-3 specialised formulas must agree with the arbitrary-`a` ones
    fn check_am3<FE, C>(g: &affine::Point<FE>)
    where
        FE: Field,
        C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveAM3,
        for<'a> &'a FE: Add<FE, Output = FE>,
        for<'a> &'a FE: Mul<FE, Output = FE>,
        for<'a> &'a FE: Sub<FE, Output = FE>,
        for<'a, 'b> &'a FE: Add<&'b FE, Output = FE>,
        for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
        for<'a, 'b> &'a FE: Sub<&'b FE, Output = FE>,
    {
        let inf = Point::<FE>::INFINITY;
        let g1 = Point::from_affine(g);
        let g2 = g1.double::<C>();

        assert_eq!(g1.double_am3::<C>().to_affine(), g2.to_affine());
        assert_eq!(
            g2.double_am3::<C>().to_affine(),
            g2.double::<C>().to_affine()
        );
        assert!(inf.double_am3::<C>().is_infinity().is_true());

        assert_eq!(
            g2.add_or_double_am3::<C>(&g1).to_affine(),
            g2.add_or_double::<C>(&g1).to_affine()
        );
        assert_eq!(g1.add_or_double_am3::<C>(&g1).to_affine(), g2.to_affine());
        assert_eq!(
            g2.add_or_double_mixed_am3::<C>(&g1).to_affine(),
            g2.add_or_double::<C>(&g1).to_affine()
        );
        assert!(g1
            .add_or_double_am3::<C>(&(-g1.clone()))
            .is_infinity()
            .is_true());

        let n: [u8; 8] = [0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xf0];
        let expected = g1.scale_ct::<C>(&n).to_affine();
        assert_eq!(g1.scale_am3_ct::<C>(&n).to_affine(), expected);
        assert_eq!(g1.scale_am3::<C>(&n).to_affine(), expected);
    }

    fn generator<FE: Field>(g: (&FE, &FE)) -> affine::Point<FE> {
        affine::Point {
            x: g.0.clone(),
            y: g.1.clone(),
        }
    }

    #[cfg(feature = "p256r1")]
    #[test]
    fn p256r1() {
        use crate::curve::sec2::p256r1::{Curve, FieldElement};
        let g = generator::<FieldElement>(Curve::generator());
        check_curve::<FieldElement, Curve>(&g);
        check_am3::<FieldElement, Curve>(&g);
    }

    #[cfg(feature = "p256k1")]
    #[test]
    fn p256k1() {
        use crate::curve::sec2::p256k1::{Curve, FieldElement};
        let g = generator::<FieldElement>(Curve::generator());
        check_curve::<FieldElement, Curve>(&g);
        check_a0::<FieldElement, Curve>(&g);
    }

    #[cfg(feature = "p384r1")]
    #[test]
    fn p384r1() {
        use crate::curve::sec2::p384r1::{Curve, FieldElement};
        let g = generator::<FieldElement>(Curve::generator());
        check_curve::<FieldElement, Curve>(&g);
        check_am3::<FieldElement, Curve>(&g);
    }
}
