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
use super::weierstrass::{WeierstrassCurve, WeierstrassCurveA0, WeierstrassCurveAM3};
use crate::mp::ct::{Choice, CtEqual, CtSelect};
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

/// This is the result return when trying to convert a projective point at infinity to its affine point
#[derive(Debug, Clone, Copy)]
pub struct AffineAtInfinity;

/// Compute the width-`w` non-adjacent form (wNAF) of the big-endian integer
/// `n`, returned as signed digits, least-significant first.
///
/// Every non-zero digit is odd and lies in `[-(2^(w-1)-1), 2^(w-1)-1]`, and no
/// two consecutive digits are both non-zero (so on average only `1/(w+1)` of
/// the digits are non-zero). This is used solely by the *variable-time* scalar
/// multiplication and is deliberately not constant-time.
fn wnaf(n: &[u8], w: u32) -> Vec<i8> {
    debug_assert!((2..=8).contains(&w));

    // Little-endian working copy of `n` with one extra byte of headroom: the
    // `k -= digit` step transiently *increases* k when the digit is negative,
    // and the spare byte guarantees the carry never runs off the end.
    let mut k = Vec::with_capacity(n.len() + 1);
    k.extend(n.iter().rev().copied());
    k.push(0);

    let width = 1i32 << w; // 2^w
    let half = 1i32 << (w - 1); // 2^(w-1)

    let mut naf = Vec::with_capacity(n.len() * 8 + 1);
    while k.iter().any(|&b| b != 0) {
        let mut digit = 0i32;
        if k[0] & 1 == 1 {
            // k mod 2^w (w <= 8, so only the low byte contributes)
            let m = (k[0] as i32) & (width - 1);
            digit = if m >= half { m - width } else { m };

            // k -= digit (signed: adds |digit| when digit < 0). The arithmetic
            // right shift turns the per-byte overflow into a sign-extended
            // carry/borrow that propagates to the next byte.
            let mut carry = -digit;
            let mut i = 0;
            while carry != 0 && i < k.len() {
                let v = k[i] as i32 + carry;
                k[i] = (v & 0xff) as u8;
                carry = v >> 8;
                i += 1;
            }
        }
        naf.push(digit as i8);

        // k >>= 1
        let mut prev = 0u8;
        for b in k.iter_mut().rev() {
            let cur = *b;
            *b = (cur >> 1) | (prev << 7);
            prev = cur & 1;
        }
    }
    naf
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
        self.z.ct_eq(&FE::ZERO)
    }
}

impl<FE> Point<FE>
where
    FE: Field,
{
    /// Returns the point at infinity
    pub const INFINITY: Self = Point {
        x: FE::ZERO,
        y: FE::ONE,
        z: FE::ZERO,
    };

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

            self.x = &self.x * &zinv;
            self.y = &self.y * &zinv;
            self.z = FE::ONE
        }
    }
}

impl<FE: Field> Point<FE> {
    pub fn add_different<'x, 'y, C: WeierstrassCurve<FieldElement = FE>>(
        &'x self,
        other: &'y Point<FE>,
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
        let z3 = C::A * &t4;
        let x3 = C::B3 * &t2;
        let z3 = &x3 + &z3;
        let x3 = &t1 - &z3;
        let z3 = &t1 + &z3;
        let y3 = &x3 * &z3;
        let t1 = t0.double();
        let t1 = t1 + &t0;
        let t2 = C::A * &t2;
        let t4 = C::B3 * &t4;
        let t1 = t1 + &t2;
        let t2 = &t0 - &t2;
        let t2 = C::A * &t2;
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
        let t2 = C::B3 * &t2;
        let z3 = &t1 + &t2;
        let t1 = t1 - &t2;
        let y3 = C::B3 * y3;
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

    pub fn add_different_am3<'x, 'y, C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveAM3>(
        &'x self,
        other: &'y Point<FE>,
    ) -> Point<FE>
    where
        for<'a> &'a FE: Add<FE, Output = FE>,
        for<'a> &'a FE: Mul<FE, Output = FE>,
        for<'a> &'a FE: Sub<FE, Output = FE>,
        for<'a, 'b> &'a FE: Add<&'b FE, Output = FE>,
        for<'a, 'b> &'a FE: Mul<&'b FE, Output = FE>,
        for<'a, 'b> &'a FE: Sub<&'b FE, Output = FE>,
    {
        // Algorithm 4 from (1) - addition formula for a=-3
        //
        // ```magma
        // ADD := function ( X1 , Y1 , Z1 , X2 , Y2 , Z2 , b )
        //     t0 := X1 * X2 ; t1 := Y1 * Y2 ; t2 := Z1 * Z2 ;
        //     t3 := X1 + Y1 ; t4 := X2 + Y2 ; t3 := t3 * t4 ;
        //     t4 := t0 + t1 ; t3 := t3 - t4 ; t4 := Y1 + Z1 ;
        //     X3 := Y2 + Z2 ; t4 := t4 * X3 ; X3 := t1 + t2 ;
        //     t4 := t4 - X3 ; X3 := X1 + Z1 ; Y3 := X2 + Z2 ;
        //     X3 := X3 * Y3 ; Y3 := t0 + t2 ; Y3 := X3 - Y3 ;
        //     Z3 := b * t2 ; X3 := Y3 - Z3 ; Z3 := X3 + X3 ;
        //     X3 := X3 + Z3 ; Z3 := t1 - X3 ; X3 := t1 + X3 ;
        //     Y3 := b * Y3 ; t1 := t2 + t2 ; t2 := t1 + t2 ;
        //     Y3 := Y3 - t2 ; Y3 := Y3 - t0 ; t1 := Y3 + Y3 ;
        //     Y3 := t1 + Y3 ; t1 := t0 + t0 ; t0 := t1 + t0 ;
        //     t0 := t0 - t2 ; t1 := t4 * Y3 ; t2 := t0 * Y3 ;
        //     Y3 := X3 * Z3 ; Y3 := Y3 + t2 ; X3 := t3 * X3 ;
        //     X3 := X3 - t1 ; Z3 := t4 * Z3 ; t1 := t3 * t0 ;
        //     Z3 := Z3 + t1 ;
        //     return X3 , Y3 , Z3 ;
        // end function ;
        // ```
        let t0 = &self.x * &other.x;
        let t1 = &self.y * &other.y;
        let t2 = &self.z * &other.z;
        let t3 = &self.x + &self.y;
        let t4 = &other.x + &other.y;
        let t3 = &t3 * &t4;
        let t4 = &t0 + &t1;
        let t3 = &t3 - &t4;
        let t4 = &self.y + &self.z;
        let x3 = &other.y + &other.z;
        let t4 = &t4 * &x3;
        let x3 = &t1 + &t2;
        let t4 = &t4 - &x3;
        let x3 = &self.x + &self.z;
        let y3 = &other.x + &other.z;
        let x3 = &x3 * &y3;
        let y3 = &t0 + &t2;
        let y3 = &x3 - &y3;
        let z3 = C::B * &t2;
        let x3 = &y3 - &z3;
        let z3 = &x3 + &x3;
        let x3 = &x3 + &z3;
        let z3 = &t1 - &x3;
        let x3 = &t1 + &x3;
        let y3 = C::B * &y3;
        let t1 = &t2 + &t2;
        let t2 = &t1 + &t2;
        let y3 = &y3 - &t2;
        let y3 = &y3 - &t0;
        let t1 = &y3 + &y3;
        let y3 = &t1 + &y3;
        let t1 = &t0 + &t0;
        let t0 = &t1 + &t0;
        let t0 = &t0 - &t2;
        let t1 = &t4 * &y3;
        let t2 = &t0 * &y3;
        let y3 = &x3 * &z3;
        let y3 = &y3 + &t2;
        let x3 = &t3 * &x3;
        let x3 = &x3 - &t1;
        let z3 = &t4 * &z3;
        let t1 = &t3 * &t0;
        let z3 = &z3 + &t1;

        Point {
            x: x3,
            y: y3,
            z: z3,
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

    /// Build the runtime fixed-base comb table from the statically embedded
    /// `sage/comb.sage` table (`src/params/comb/<curve>.rs`).
    ///
    /// `table[i]` holds the 15 affine points `1·B, 2·B, …, 15·B` (where
    /// `B = 16^i·G` for window `i`) as `(x, y)` coordinate byte arrays. The
    /// returned per-window arrays place those at indices `1..=15`, with index
    /// `0` set to the point at infinity so a zero digit selects the neutral
    /// element. `parse` is the coordinate decoder (`FieldElement::from_bytes`).
    ///
    /// The result is heap-allocated (`Box`): the table is large for the bigger
    /// curves (e.g. p521 needs `132·16` points ≈ 450 KiB), so building it on the
    /// stack would overflow it in unoptimized builds. We grow a `Vec` one window
    /// at a time — only a single 16-point window is ever on the stack — and then
    /// reuse that allocation as the boxed array (no copy of the bulk data).
    pub(crate) fn build_comb_table<const NW: usize, const FS: usize>(
        table: &[[([u8; FS], [u8; FS]); 15]; NW],
        parse: fn(&[u8; FS]) -> FE,
    ) -> Box<[[Point<FE>; 16]; NW]> {
        let mut windows: Vec<[Point<FE>; 16]> = Vec::with_capacity(NW);
        for w in 0..NW {
            let mut window: [Point<FE>; 16] = core::array::from_fn(|_| Self::INFINITY);
            for (slot, (x, y)) in window.iter_mut().skip(1).zip(table[w].iter()) {
                *slot = Point {
                    x: parse(x),
                    y: parse(y),
                    z: FE::ONE,
                };
            }
            windows.push(window);
        }

        // conversion should not reallocate since we pinned the capacity
        <Box<[[Point<FE>; 16]; NW]>>::try_from(windows.into_boxed_slice())
            .ok()
            .expect("comb table window count matches NW")
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
    pub fn double<C: WeierstrassCurve<FieldElement = FE>>(&self) -> Self {
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
        let x3 = C::A * &z3;
        let y3 = C::B3 * &t2;
        let y3 = &x3 + &y3;
        let x3 = &t1 - &y3;
        let y3 = &t1 + &y3;
        let y3 = &x3 * &y3;
        let x3 = &t3 * &x3;
        let z3 = C::B3 * &z3;
        let t2 = C::A * &t2;
        let t3 = &t0 - &t2;
        let t3 = C::A * &t3;
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
    pub fn double_a0<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveA0>(&self) -> Self {
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
        let t2 = C::B3 * &t2;
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

    #[inline]
    pub fn double_am3<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveAM3>(&self) -> Self {
        // Algorithm 6 from (1) - doubling formula for a=-3
        //
        // ```magma
        // DBL := function (X ,Y ,Z , b )
        //     t0 := X ^2; t1 := Y ^2; t2 := Z ^2;
        //     t3 := X * Y ; t3 := t3 + t3 ; Z3 := X * Z ;
        //     Z3 := Z3 + Z3 ; Y3 := b * t2 ; Y3 := Y3 - Z3 ;
        //     X3 := Y3 + Y3 ; Y3 := X3 + Y3 ; X3 := t1 - Y3 ;
        //     Y3 := t1 + Y3 ; Y3 := X3 * Y3 ; X3 := X3 * t3 ;
        //     t3 := t2 + t2 ; t2 := t2 + t3 ; Z3 := b * Z3 ;
        //     Z3 := Z3 - t2 ; Z3 := Z3 - t0 ; t3 := Z3 + Z3 ;
        //     Z3 := Z3 + t3 ; t3 := t0 + t0 ; t0 := t3 + t0 ;
        //     t0 := t0 - t2 ; t0 := t0 * Z3 ; Y3 := Y3 + t0 ;
        //     t0 := Y * Z ; t0 := t0 + t0 ; Z3 := t0 * Z3 ;
        //     X3 := X3 - Z3 ; Z3 := t0 * t1 ; Z3 := Z3 + Z3 ;
        //     Z3 := Z3 + Z3 ;
        //     return X3 , Y3 , Z3 ;
        // end function ;
        // ```
        let t0 = self.x.square();
        let t1 = self.y.square();
        let t2 = self.z.square();
        let t3 = &self.x * &self.y;
        let t3 = &t3 + &t3;
        let z3 = &self.x * &self.z;
        let z3 = &z3 + &z3;
        let y3 = C::B * &t2;
        let y3 = &y3 - &z3;
        let x3 = &y3 + &y3;
        let y3 = &x3 + &y3;
        let x3 = &t1 - &y3;
        let y3 = &t1 + &y3;
        let y3 = &x3 * &y3;
        let x3 = &x3 * &t3;
        let t3 = &t2 + &t2;
        let t2 = &t2 + &t3;
        let z3 = C::B * &z3;
        let z3 = &z3 - &t2;
        let z3 = &z3 - &t0;
        let t3 = &z3 + &z3;
        let z3 = &z3 + &t3;
        let t3 = &t0 + &t0;
        let t0 = &t3 + &t0;
        let t0 = &t0 - &t2;
        let t0 = &t0 * &z3;
        let y3 = &y3 + &t0;
        let t0 = &self.y * &self.z;
        let t0 = &t0 + &t0;
        let z3 = &t0 * &z3;
        let x3 = &x3 - &z3;
        let z3 = &t0 * &t1;
        let z3 = &z3 + &z3;
        let z3 = &z3 + &z3;

        Point {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    pub fn to_affine(&self) -> Option<affine::Point<FE>> {
        if self.z == FE::ONE {
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

    /// Width of the signed window used by the variable-time wNAF scalar
    /// multiplication. `W = 5` keeps the odd-multiple table small (`2^(W-2)`
    /// points) while giving an average non-zero digit density of `1/(W+1)`,
    /// which is close to optimal for ~192–521 bit scalars.
    const WNAF_W: u32 = 5;

    /// Additive inverse of the point: `-(X:Y:Z) = (X:-Y:Z)`.
    #[inline]
    fn negate(&self) -> Self {
        Point {
            x: self.x.clone(),
            y: -self.y.clone(),
            z: self.z.clone(),
        }
    }

    /// Variable-time scalar multiplication `n * self` using a width-`W` wNAF.
    ///
    /// Compared with a plain double-and-add this keeps the same number of
    /// doublings but cuts the additions from ~`nbits/2` down to ~`nbits/(W+1)`
    /// by using a signed sparse recoding of the scalar. The running time
    /// depends on the scalar, so this must only be used when `n` is public.
    ///
    /// `add_different` is the *complete* addition formula, so the `q == self`
    /// and `q == infinity` cases need no special handling.
    #[inline]
    fn scalar_mul_wnaf<C: WeierstrassCurve<FieldElement = FE>>(&self, n: &[u8]) -> Self {
        let naf = wnaf(n, Self::WNAF_W);
        // table[i] = (2i+1) * self, i.e. the odd multiples 1·P, 3·P, … of self
        let dbl = self.double::<C>();
        let tlen = 1usize << (Self::WNAF_W - 2);
        let mut table: Vec<Point<FE>> = Vec::with_capacity(tlen);
        table.push(self.clone());
        for i in 1..tlen {
            table.push(table[i - 1].add_different::<C>(&dbl));
        }

        let mut q = Self::INFINITY;
        for &d in naf.iter().rev() {
            q = q.double::<C>();
            if d > 0 {
                q = q.add_different::<C>(&table[(d as usize) >> 1]);
            } else if d < 0 {
                q = q.add_different::<C>(&table[((-(d as i32)) as usize) >> 1].negate());
            }
        }
        q
    }

    #[inline]
    fn scalar_mul_wnaf_a0<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveA0>(
        &self,
        n: &[u8],
    ) -> Self {
        let naf = wnaf(n, Self::WNAF_W);
        let dbl = self.double_a0::<C>();
        let tlen = 1usize << (Self::WNAF_W - 2);
        let mut table: Vec<Point<FE>> = Vec::with_capacity(tlen);
        table.push(self.clone());
        for i in 1..tlen {
            table.push(table[i - 1].add_different_a0::<C>(&dbl));
        }

        let mut q = Self::INFINITY;
        for &d in naf.iter().rev() {
            q = q.double_a0::<C>();
            if d > 0 {
                q = q.add_different_a0::<C>(&table[(d as usize) >> 1]);
            } else if d < 0 {
                q = q.add_different_a0::<C>(&table[((-(d as i32)) as usize) >> 1].negate());
            }
        }
        q
    }

    #[inline]
    fn scalar_mul_wnaf_am3<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveAM3>(
        &self,
        n: &[u8],
    ) -> Self {
        let naf = wnaf(n, Self::WNAF_W);
        let dbl = self.double_am3::<C>();
        let tlen = 1usize << (Self::WNAF_W - 2);
        let mut table: Vec<Point<FE>> = Vec::with_capacity(tlen);
        table.push(self.clone());
        for i in 1..tlen {
            table.push(table[i - 1].add_different_am3::<C>(&dbl));
        }

        let mut q = Self::INFINITY;
        for &d in naf.iter().rev() {
            q = q.double_am3::<C>();
            if d > 0 {
                q = q.add_different_am3::<C>(&table[(d as usize) >> 1]);
            } else if d < 0 {
                q = q.add_different_am3::<C>(&table[((-(d as i32)) as usize) >> 1].negate());
            }
        }
        q
    }

    pub fn scale<C: WeierstrassCurve<FieldElement = FE>>(&self, n: &[u8]) -> Self {
        self.scalar_mul_wnaf::<C>(n)
    }

    pub fn scale_a0<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveA0>(
        &self,
        n: &[u8],
    ) -> Self {
        self.scalar_mul_wnaf_a0::<C>(n)
    }

    pub fn scale_am3<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveAM3>(
        &self,
        n: &[u8],
    ) -> Self {
        self.scalar_mul_wnaf_am3::<C>(n)
    }

    /// Constant-time scalar multiplication using a fixed 4-bit window.
    ///
    /// The window value is read from a precomputed table with a constant-time
    /// scan, and the addition relies on the complete formula so a zero window
    /// needs no special case.
    fn scalar_mul_fixed_window<C: WeierstrassCurve<FieldElement = FE>>(&self, n: &[u8]) -> Self {
        // table[d] = d * self, for d in 0..16
        let mut table: [Point<FE>; 16] = core::array::from_fn(|_| Self::INFINITY);
        table[1] = self.clone();
        table[2] = self.double::<C>();
        for d in 3..16 {
            let next = table[d - 1].add_different::<C>(self);
            table[d] = next;
        }

        let mut q = Self::INFINITY;
        for byte in n.iter() {
            for &index in &[byte >> 4, byte & 0x0f] {
                q = q.double::<C>().double::<C>().double::<C>().double::<C>();
                let selected = Self::select_from_table(&table, index);
                q = q.add_different::<C>(&selected);
            }
        }
        q
    }

    /// Constant-time fixed 4-bit window scalar multiplication for a=0 curves.
    /// See [`Self::scalar_mul_fixed_window`].
    fn scalar_mul_fixed_window_a0<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveA0>(
        &self,
        n: &[u8],
    ) -> Self {
        let mut table: [Point<FE>; 16] = core::array::from_fn(|_| Self::INFINITY);
        table[1] = self.clone();
        table[2] = self.double_a0::<C>();
        for d in 3..16 {
            let next = table[d - 1].add_different_a0::<C>(self);
            table[d] = next;
        }

        let mut q = Self::INFINITY;
        for byte in n.iter() {
            for &index in &[byte >> 4, byte & 0x0f] {
                q = q
                    .double_a0::<C>()
                    .double_a0::<C>()
                    .double_a0::<C>()
                    .double_a0::<C>();
                let selected = Self::select_from_table(&table, index);
                q = q.add_different_a0::<C>(&selected);
            }
        }
        q
    }

    /// Constant-time fixed 4-bit window scalar multiplication for a=-3 curves.
    /// See [`Self::scalar_mul_fixed_window`].
    fn scalar_mul_fixed_window_am3<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveAM3>(
        &self,
        n: &[u8],
    ) -> Self {
        let mut table: [Point<FE>; 16] = core::array::from_fn(|_| Self::INFINITY);
        table[1] = self.clone();
        table[2] = self.double_am3::<C>();
        for d in 3..16 {
            let next = table[d - 1].add_different_am3::<C>(self);
            table[d] = next;
        }

        let mut q = Self::INFINITY;
        for byte in n.iter() {
            for &index in &[byte >> 4, byte & 0x0f] {
                q = q
                    .double_am3::<C>()
                    .double_am3::<C>()
                    .double_am3::<C>()
                    .double_am3::<C>();
                let selected = Self::select_from_table(&table, index);
                q = q.add_different_am3::<C>(&selected);
            }
        }
        q
    }

    /// Constant-time scalar multiplication (default). See
    /// [`Self::scalar_mul_fixed_window`].
    pub fn scale_ct<C: WeierstrassCurve<FieldElement = FE>>(&self, n: &[u8]) -> Self {
        self.scalar_mul_fixed_window::<C>(n)
    }

    /// Constant-time scalar multiplication for a=0 curves (default).
    pub fn scale_a0_ct<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveA0>(
        &self,
        n: &[u8],
    ) -> Self {
        self.scalar_mul_fixed_window_a0::<C>(n)
    }

    /// Constant-time scalar multiplication for a=-3 curves (default).
    pub fn scale_am3_ct<C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveAM3>(
        &self,
        n: &[u8],
    ) -> Self {
        self.scalar_mul_fixed_window_am3::<C>(n)
    }

    /// Constant-time fixed-base scalar multiplication from a precomputed comb
    /// table (see [`Self::build_comb_table`]): computes `n · G`, where `G` is the
    /// generator the table was built from.
    ///
    /// The table already folds in the per-window `16^i` weights, so this needs
    /// no point doublings: just one constant-time table lookup and one complete
    /// addition per 4-bit window of the scalar. The digit only ever feeds the
    /// constant-time table scan, so the running time is independent of `n`.
    pub fn mul_base_table<const NW: usize, C: WeierstrassCurve<FieldElement = FE>>(
        tables: &[[Self; 16]; NW],
        n: &[u8],
    ) -> Self {
        assert_eq!(tables.len(), n.len() * 2, "comb table size mismatch");
        let mut q = Self::INFINITY;
        for (i, window) in tables.iter().enumerate() {
            let byte = n[n.len() - 1 - (i / 2)];
            let digit = if i % 2 == 0 { byte & 0x0f } else { byte >> 4 };
            let selected = Self::select_from_table(window, digit);
            q = q.add_different::<C>(&selected);
        }
        q
    }

    /// Fixed-base comb scalar multiplication for a=0 curves. See
    /// [`Self::mul_base_table`].
    pub fn mul_base_table_a0<
        const NW: usize,
        C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveA0,
    >(
        tables: &[[Point<FE>; 16]; NW],
        n: &[u8],
    ) -> Self {
        assert_eq!(tables.len(), n.len() * 2, "comb table size mismatch");
        let mut q = Self::INFINITY;
        for (i, window) in tables.iter().enumerate() {
            let byte = n[n.len() - 1 - (i / 2)];
            let digit = if i % 2 == 0 { byte & 0x0f } else { byte >> 4 };
            let selected = Self::select_from_table(window, digit);
            q = q.add_different_a0::<C>(&selected);
        }
        q
    }

    /// Fixed-base comb scalar multiplication for a=-3 curves. See
    /// [`Self::mul_base_table`].
    pub fn mul_base_table_am3<
        const NW: usize,
        C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveAM3,
    >(
        tables: &[[Point<FE>; 16]; NW],
        n: &[u8],
    ) -> Self {
        assert_eq!(tables.len(), n.len() * 2, "comb table size mismatch");
        let mut q = Self::INFINITY;
        for (i, window) in tables.iter().enumerate() {
            let byte = n[n.len() - 1 - (i / 2)];
            let digit = if i % 2 == 0 { byte & 0x0f } else { byte >> 4 };
            let selected = Self::select_from_table(window, digit);
            q = q.add_different_am3::<C>(&selected);
        }
        q
    }

    /// Add two points, correctly handling every case (`self == other`,
    /// `self == -other` and the point at infinity).
    ///
    /// This relies on the *complete* addition formula (Algorithm 1), which is
    /// exception-free for prime order curves, so unlike a naive implementation
    /// it does not need to compare the two points and branch to a dedicated
    /// doubling routine.
    #[inline]
    pub fn add_or_double<'b, C: WeierstrassCurve<FieldElement = FE>>(
        &self,
        other: &'b Point<FE>,
    ) -> Point<FE> {
        self.add_different::<C>(other)
    }

    #[inline]
    pub fn add_or_double_a0<'b, C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveA0>(
        &self,
        other: &'b Point<FE>,
    ) -> Point<FE> {
        self.add_different_a0::<C>(other)
    }

    #[inline]
    pub fn add_or_double_am3<'b, C: WeierstrassCurve<FieldElement = FE> + WeierstrassCurveAM3>(
        &self,
        other: &'b Point<FE>,
    ) -> Point<FE> {
        self.add_different_am3::<C>(other)
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
