//! Jubjub
//!
//! Jubjub is the twisted Edwards curve
//!
//! ```text
//! -x^2 + y^2 = 1 + d*x^2*y^2  ,  d = -(10240/10241)
//! ```
//!
//! defined over the BLS12-381 scalar field `Fr` (i.e. its coordinates are
//! [`bls12_381::Scalar`](crate::curve::bls12_381::Scalar) elements, re-exported
//! here as [`FieldElement`]). It is a complete `a = -1` twisted Edwards curve,
//! so the group law has no exceptional cases.
//!
//! The full curve group has order `8 * r`; [`Point`] implements the prime-order
//! subgroup of 252-bit order `r`, whose scalars form the field
//! [`Scalar`](scalar::Scalar). This mirrors the edwards25519 implementation in
//! [`crate::curve::curve25519`], reusing the same complete `a = -1` addition and
//! doubling formulas.

use crate::curve::edwards::{EdwardsCurve, EdwardsCurveAM1};
use crate::curve::field::Sign;
use crate::mp::ct::{Choice, CtEqual, CtSelect};
use crate::params::jubjub::{A_BYTES, D2_BYTES, D_BYTES, GT_BYTES, GX_BYTES, GY_BYTES};
#[cfg(feature = "table")]
use crate::params::jubjub::{COMB_TABLE, COMB_WINDOWS};
#[cfg(feature = "table")]
use core::convert::TryFrom;
use core::ops::{Add, Mul, Neg, Sub};

pub mod scalar;

pub use scalar::Scalar;

/// The Jubjub base field: the BLS12-381 scalar field `Fr`, over which the curve
/// coordinates live.
pub use crate::curve::bls12_381::Scalar as FieldElement;

/// The twisted Edwards curve `-x^2 + y^2 = 1 + d*x^2*y^2` (Jubjub).
#[derive(Debug, Clone, Copy)]
pub struct EdCurve;

impl EdwardsCurve for EdCurve {
    type FieldElement = FieldElement;
    const A: FieldElement = FieldElement::from_bytes_unchecked_be(&A_BYTES);
    const D: FieldElement = FieldElement::from_bytes_unchecked_be(&D_BYTES);
    const D2: FieldElement = FieldElement::from_bytes_unchecked_be(&D2_BYTES);
}
impl EdwardsCurveAM1 for EdCurve {}

// ***********************************************************************
// Twisted Edwards point (prime-order subgroup)
// ***********************************************************************

/// Point on the Jubjub twisted Edwards curve, in extended homogeneous
/// coordinates `(X:Y:Z:T)` with `x = X/Z`, `y = Y/Z` and `T = XY/Z`.
///
/// Uses the complete addition formulas valid for `a = -1`
/// (see [`EdwardsCurveAM1`]), so there are no exceptional cases.
#[derive(Clone, Debug)]
pub struct Point {
    x: FieldElement,
    y: FieldElement,
    z: FieldElement,
    t: FieldElement,
}

impl Point {
    /// Neutral element of the group (the affine point `(0, 1)`)
    pub const IDENTITY: Self = Point {
        x: FieldElement::zero(),
        y: FieldElement::one(),
        z: FieldElement::one(),
        t: FieldElement::zero(),
    };

    /// Curve generator (base point) of the prime-order subgroup
    pub const GENERATOR: Self = Point {
        x: FieldElement::from_bytes_unchecked_be(&GX_BYTES),
        y: FieldElement::from_bytes_unchecked_be(&GY_BYTES),
        z: FieldElement::one(),
        t: FieldElement::from_bytes_unchecked_be(&GT_BYTES),
    };

    fn from_affine(x: &FieldElement, y: &FieldElement) -> Self {
        Point {
            x: x.clone(),
            y: y.clone(),
            z: FieldElement::one(),
            t: x * y,
        }
    }

    /// Try to create a point from its affine coordinates `(x, y)`, returning
    /// `None` if the point is not on the curve.
    pub fn from_coordinate(x: &FieldElement, y: &FieldElement) -> Option<Self> {
        // a*x^2 + y^2 == 1 + d*x^2*y^2 , with a = -1
        let xx = x.square();
        let yy = y.square();
        let lhs = &yy - &xx;
        let rhs = &FieldElement::one() + &(&EdCurve::D * &(&xx * &yy));
        if lhs.ct_eq(&rhs).is_true() {
            Some(Point::from_affine(x, y))
        } else {
            None
        }
    }

    /// Return the affine coordinates `(x, y)` of the point
    pub fn to_affine(&self) -> (FieldElement, FieldElement) {
        let zinv = self.z.inverse();
        (&self.x * &zinv, &self.y * &zinv)
    }

    /// Point doubling using the dedicated `a = -1` formula.
    pub fn double(&self) -> Point {
        let a = self.x.square();
        let b = self.y.square();
        let c = self.z.square().double(); // 2*Z^2
        let d = -&a; // a*A with a = -1
        let xy = &self.x + &self.y;
        let e = &xy.square() - &(&a + &b); // (X+Y)^2 - A - B
        let g = &d + &b;
        let f = &g - &c;
        let h = &d - &b;
        Point {
            x: &e * &f,
            y: &g * &h,
            z: &f * &g,
            t: &e * &h,
        }
    }

    /// Complete point addition (valid for every pair of points since `a = -1`).
    pub fn add(&self, other: &Point) -> Point {
        let aa = &(&self.y - &self.x) * &(&other.y - &other.x);
        let bb = &(&self.y + &self.x) * &(&other.y + &other.x);
        let cc = &(&EdCurve::D2 * &self.t) * &other.t; // 2*d*T1*T2
        let dd = (&self.z * &other.z).double(); // 2*Z1*Z2
        let e = &bb - &aa;
        let f = &dd - &cc;
        let g = &dd + &cc;
        let h = &bb + &aa;
        Point {
            x: &e * &f,
            y: &g * &h,
            z: &f * &g,
            t: &e * &h,
        }
    }

    fn negate(&self) -> Point {
        Point {
            x: -&self.x,
            y: self.y.clone(),
            z: self.z.clone(),
            t: -&self.t,
        }
    }

    /// Constant-time scalar multiplication by a big-endian scalar.
    ///
    /// A simple double-and-add over the complete addition formula: the optional
    /// addition is selected in constant time so the running time does not depend
    /// on the scalar bits.
    pub fn scale_bytes(&self, k_be: &[u8]) -> Point {
        let mut q = Point::IDENTITY;
        for &byte in k_be.iter() {
            for i in (0..8).rev() {
                q = q.double();
                let bit = (byte >> i) & 1;
                let added = Point::add(&q, self);
                q = Point::ct_select(Choice((bit as u64) & 1), &added, &q);
            }
        }
        q
    }

    /// Constant-time scalar multiplication by a [`Scalar`].
    pub fn scale(&self, k: &Scalar) -> Point {
        self.scale_bytes(&k.to_bytes_be())
    }

    /// Point compression: the y-coordinate together with the sign of x.
    pub fn compress(&self) -> (FieldElement, Sign) {
        let (x, y) = self.to_affine();
        (y, x.sign())
    }

    /// Point decompression: recover the point from a y-coordinate and the sign
    /// of x. Returns `None` if no curve point has this y-coordinate.
    ///
    /// Note that the base-field square root ([`FieldElement::sqrt`]) is not
    /// constant-time, so decompression should only be done on public values.
    pub fn decompress(y: &FieldElement, x_sign: Sign) -> Option<Self> {
        let one = FieldElement::one();
        let yy = y.square();
        // x^2 = (y^2 - 1) / (d*y^2 + 1)
        let u = &yy - &one;
        let v = &(&EdCurve::D * &yy) + &one;
        if v.is_zero() {
            return None;
        }
        let x2 = &u * &v.inverse();
        match x2.sqrt().into_option() {
            None => None,
            Some(x) => {
                let x = if x.sign() == x_sign { x } else { -x };
                Some(Point::from_affine(&x, y))
            }
        }
    }

    /// Constant-time fixed-base scalar multiplication: `scalar · B`, where `B`
    /// is the curve generator.
    ///
    /// With the `table` feature this uses a precomputed comb table for the
    /// generator (built once on first use), which is several times faster than
    /// the general [`Point::scale`] since it needs no point doublings: just one
    /// constant-time table lookup and one complete addition per 4-bit window of
    /// the scalar. Without the feature it falls back to `scale`.
    #[cfg(feature = "table")]
    pub fn mul_base(scalar: &Scalar) -> Point {
        let tables = generator_comb();
        let n = scalar.to_bytes_be(); // big-endian: 32 bytes = 64 nibbles
        let mut q = Point::IDENTITY;
        for (i, window) in tables.iter().enumerate() {
            let byte = n[n.len() - 1 - (i / 2)];
            let digit = if i % 2 == 0 { byte & 0x0f } else { byte >> 4 };
            let selected = Self::select_from_table(window, digit);
            q = q.add(&selected);
        }
        q
    }

    /// Fixed-base scalar multiplication `scalar · B` (no precomputation).
    #[cfg(not(feature = "table"))]
    pub fn mul_base(scalar: &Scalar) -> Point {
        Point::GENERATOR.scale(scalar)
    }

    /// Constant-time lookup of `table[index]`: the whole window is scanned so
    /// the memory access pattern does not depend on the (secret) `index`.
    #[cfg(feature = "table")]
    fn select_from_table(table: &[Point; 16], index: u8) -> Point {
        let mut acc = Point::IDENTITY;
        for (j, t) in table.iter().enumerate() {
            let take = (j as u64).ct_eq(&(index as u64));
            acc = Point::ct_select(take, t, &acc);
        }
        acc
    }
}

/// The fixed-base comb table for the generator, parsed once on first use from
/// the statically embedded [`COMB_TABLE`] (generated by `sage/comb.sage`).
///
/// `COMB_TABLE[i][d]` holds the affine `(x, y)` of `(d + 1) · 16^i · B`; the
/// runtime table places those at window indices `1..=15`, with index `0` set to
/// the identity so a zero digit selects the neutral element. `n · B` is then one
/// constant-time table lookup and one complete addition per 4-bit window of `n`,
/// with no point doublings.
#[cfg(feature = "table")]
fn generator_comb() -> &'static [[Point; 16]; COMB_WINDOWS] {
    static V: std::sync::OnceLock<alloc::boxed::Box<[[Point; 16]; COMB_WINDOWS]>> =
        std::sync::OnceLock::new();
    &**V.get_or_init(build_comb_table)
}

#[cfg(feature = "table")]
fn build_comb_table() -> alloc::boxed::Box<[[Point; 16]; COMB_WINDOWS]> {
    let mut windows: alloc::vec::Vec<[Point; 16]> = alloc::vec::Vec::with_capacity(COMB_WINDOWS);
    for row in COMB_TABLE.iter() {
        let mut window: [Point; 16] = core::array::from_fn(|_| Point::IDENTITY);
        for (slot, (x, y)) in window.iter_mut().skip(1).zip(row.iter()) {
            *slot = Point::from_affine(
                &FieldElement::from_bytes_unchecked_be(x),
                &FieldElement::from_bytes_unchecked_be(y),
            );
        }
        windows.push(window);
    }
    <alloc::boxed::Box<[[Point; 16]; COMB_WINDOWS]>>::try_from(windows.into_boxed_slice())
        .ok()
        .expect("comb window count matches COMB_WINDOWS")
}

impl CtSelect for Point {
    fn ct_select(cond: Choice, a: &Point, b: &Point) -> Point {
        Point {
            x: FieldElement::ct_select(cond, &a.x, &b.x),
            y: FieldElement::ct_select(cond, &a.y, &b.y),
            z: FieldElement::ct_select(cond, &a.z, &b.z),
            t: FieldElement::ct_select(cond, &a.t, &b.t),
        }
    }
    fn ct_assign(&mut self, cond: Choice, other: &Point) {
        *self = Self::ct_select(cond, other, self);
    }
}

impl PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        // (X1:Y1:Z1) == (X2:Y2:Z2) iff X1 * Z2 == X2 * Z1 and Y1 * Z2 == Y2 * Z1
        let x1z2 = &self.x * &other.z;
        let x2z1 = &other.x * &self.z;
        let y1z2 = &self.y * &other.z;
        let y2z1 = &other.y * &self.z;
        (x1z2.ct_eq(&x2z1) & y1z2.ct_eq(&y2z1)).is_true()
    }
}
impl Eq for Point {}

impl Neg for Point {
    type Output = Point;
    fn neg(self) -> Point {
        self.negate()
    }
}
impl<'a> Neg for &'a Point {
    type Output = Point;
    fn neg(self) -> Point {
        self.negate()
    }
}

impl<'a, 'b> Add<&'b Point> for &'a Point {
    type Output = Point;
    fn add(self, other: &'b Point) -> Point {
        Point::add(self, other)
    }
}
impl<'b> Add<&'b Point> for Point {
    type Output = Point;
    fn add(self, other: &'b Point) -> Point {
        &self + other
    }
}
impl<'a> Add<Point> for &'a Point {
    type Output = Point;
    fn add(self, other: Point) -> Point {
        self + &other
    }
}
impl Add<Point> for Point {
    type Output = Point;
    fn add(self, other: Point) -> Point {
        &self + &other
    }
}
impl<'a, 'b> Sub<&'b Point> for &'a Point {
    type Output = Point;
    fn sub(self, other: &'b Point) -> Point {
        Point::add(self, &other.negate())
    }
}
impl<'b> Sub<&'b Point> for Point {
    type Output = Point;
    fn sub(self, other: &'b Point) -> Point {
        &self - other
    }
}
impl<'a> Sub<Point> for &'a Point {
    type Output = Point;
    fn sub(self, other: Point) -> Point {
        self - &other
    }
}
impl Sub<Point> for Point {
    type Output = Point;
    fn sub(self, other: Point) -> Point {
        &self - &other
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a Point {
    type Output = Point;
    fn mul(self, k: &'b Scalar) -> Point {
        self.scale(k)
    }
}
impl<'b> Mul<&'b Scalar> for Point {
    type Output = Point;
    fn mul(self, k: &'b Scalar) -> Point {
        self.scale(k)
    }
}
impl<'a, 'b> Mul<&'b Point> for &'a Scalar {
    type Output = Point;
    fn mul(self, p: &'b Point) -> Point {
        p.scale(self)
    }
}

impl crate::curve::group::CurveGroup for Point {
    type Scalar = Scalar;

    const IDENTITY: Self = Point::IDENTITY;
    const GENERATOR: Self = Point::GENERATOR;

    fn double(&self) -> Self {
        Point::double(self)
    }
    fn mul_base(scalar: &Scalar) -> Self {
        Point::mul_base(scalar)
    }
    fn mul_vartime(&self, scalar: &Scalar) -> Self {
        self.scale(scalar)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::jubjub::ORDER_BYTES;

    // in-module test vectors are big-endian
    fn fe(bytes: &[u8; 32]) -> FieldElement {
        FieldElement::from_bytes_be(bytes).expect("valid field element")
    }

    #[test]
    fn generator_on_curve() {
        let (x, y) = Point::GENERATOR.to_affine();
        assert!(Point::from_coordinate(&x, &y).is_some());
        // the stored extended T must equal x*y
        assert_eq!(&x * &y, fe(&GT_BYTES));
    }

    #[test]
    fn identity_laws() {
        let g = Point::GENERATOR;
        assert_eq!(&g + &Point::IDENTITY, g);
        assert_eq!(&Point::IDENTITY + &g, g);
        assert_eq!(&g + &(-&g), Point::IDENTITY);
    }

    #[test]
    fn double_matches_add() {
        let g = Point::GENERATOR;
        assert_eq!(g.double(), &g + &g);
        assert_eq!(g.scale_bytes(&[2]), g.double());
        assert_eq!(g.scale_bytes(&[3]), &g.double() + &g);
    }

    #[test]
    fn generator_has_order_r() {
        // r * G == identity, so the generator lies in the prime-order subgroup
        assert_eq!(Point::GENERATOR.scale_bytes(&ORDER_BYTES), Point::IDENTITY);
    }

    #[test]
    fn scalar_linearity() {
        let g = Point::GENERATOR;
        for k in [1u8, 2, 5, 9, 200] {
            let lhs = g.scale_bytes(&[k]);
            let mut rhs = Point::IDENTITY;
            for _ in 0..k {
                rhs = &rhs + &g;
            }
            assert_eq!(lhs, rhs, "k={}", k);
        }
    }

    #[test]
    fn compress_roundtrip() {
        let g = Point::GENERATOR;
        for k in [1u8, 2, 3, 7, 50] {
            let p = g.scale_bytes(&[k]);
            let (y, sign) = p.compress();
            let q = Point::decompress(&y, sign).expect("decompress");
            assert_eq!(p, q, "k={}", k);
        }
    }

    #[test]
    fn mul_base_matches_scale() {
        assert_eq!(Point::mul_base(&Scalar::from_u64(0)), Point::IDENTITY);
        for k in [1u64, 2, 9, 16, 255, 256, 1234567, 0x0123_4567_89ab_cdef] {
            let s = Scalar::from_u64(k);
            assert_eq!(Point::mul_base(&s), Point::GENERATOR.scale(&s), "k={}", k);
        }
        let mut s = Scalar::from_u64(0x0123_4567_89ab_cdef);
        for _ in 0..5 {
            s = s.square();
        }
        assert_eq!(Point::mul_base(&s), Point::GENERATOR.scale(&s));
    }

    // independently (sage) computed affine coordinates of k*G, big-endian.
    #[rustfmt::skip]
    const KAT: &[(u64, [u8; 32], [u8; 32])] = &[
        (2, [
            0x5f,0x79,0x3e,0xc6,0x3b,0x98,0xd1,0x9e,0xde,0xc4,0x26,0x04,0x1a,0xc6,0x30,0x39,
            0xf5,0x69,0x29,0xd3,0xa1,0x41,0xaf,0xff,0x61,0xd9,0xa5,0x87,0x5a,0xd0,0x6a,0xb5,
        ], [
            0x51,0x6d,0xef,0x4d,0x9d,0x10,0x0e,0xab,0x4c,0x21,0xe3,0x56,0xd2,0x0c,0x48,0x51,
            0x07,0xaf,0xee,0x33,0xf5,0xde,0xae,0x32,0x53,0xb4,0xf1,0x00,0x58,0xc0,0xdf,0x97,
        ]),
        (3, [
            0x5d,0xfd,0xd6,0x3e,0xc4,0xbf,0xf0,0x0e,0x6c,0xa7,0x57,0xed,0x97,0x41,0x7d,0x04,
            0xaa,0xbe,0x4e,0x8f,0xf1,0x6d,0x22,0x73,0xca,0xab,0xd2,0x83,0x5b,0xa5,0x55,0x30,
        ], [
            0x4f,0x42,0x81,0xba,0x2e,0x87,0x79,0x06,0x0a,0x7a,0x20,0xcf,0x6f,0x5e,0xb5,0xa3,
            0xd3,0xbf,0x8e,0x9f,0xbd,0xa4,0xcd,0xe2,0x36,0x70,0xa5,0x12,0xcd,0x0f,0x31,0x5a,
        ]),
        (5, [
            0x4e,0xca,0x40,0xcb,0x8c,0x18,0x72,0xe9,0x24,0xda,0x41,0x89,0x83,0x71,0x00,0xf4,
            0x06,0x03,0x60,0xeb,0x56,0xb0,0xe1,0x32,0xca,0xb0,0xfe,0x51,0x77,0x5c,0x2f,0x09,
        ], [
            0x42,0x7f,0x08,0x31,0xfa,0x70,0xda,0x08,0x65,0xa8,0x52,0x57,0xda,0xd6,0x65,0x84,
            0x3f,0x41,0x81,0xb2,0xb0,0xca,0xe3,0xee,0x5c,0x35,0xe3,0xfe,0x7a,0xa8,0xda,0x82,
        ]),
        (100, [
            0x08,0x9c,0x7f,0xf8,0x59,0xae,0x59,0xf7,0xcf,0x63,0x63,0x21,0x43,0x20,0x08,0xc6,
            0x47,0x71,0x4f,0x64,0xc7,0x0e,0x3a,0x9a,0x61,0xd5,0xc5,0x5d,0x0f,0xb8,0xc8,0xa5,
        ], [
            0x6f,0xcc,0x48,0xdb,0xab,0x1e,0x8f,0x9a,0x38,0x66,0x4d,0x8a,0xee,0x33,0x5e,0x49,
            0xaf,0x02,0x9f,0xc3,0xa1,0xec,0x27,0xd2,0x91,0xc6,0x72,0xe0,0xd7,0x24,0x88,0x10,
        ]),
        (1000, [
            0x1a,0xde,0x07,0x63,0x99,0xac,0xd2,0x7c,0x81,0xb7,0x4e,0xca,0x75,0x10,0x56,0x47,
            0x48,0x35,0x0a,0x11,0xc7,0xed,0x99,0x52,0xa4,0x47,0x19,0xb1,0xd1,0x8b,0x8e,0x4b,
        ], [
            0x07,0xef,0x72,0x12,0xac,0x81,0x22,0xed,0xa0,0x5f,0x6a,0xb9,0xe2,0xdd,0x88,0xcd,
            0xa3,0x11,0xda,0x6d,0xcf,0x15,0xce,0x7c,0x41,0x29,0x65,0x1c,0x4b,0xb0,0xd8,0x76,
        ]),
        (0x0123_4567_89ab_cdef, [
            0x23,0xdf,0x23,0x46,0x42,0x7b,0x16,0x18,0xf8,0xfe,0x65,0xdb,0x9f,0xa1,0xcc,0x82,
            0x26,0xf0,0xb9,0x79,0xc5,0xc5,0x9e,0xfd,0xc7,0xdf,0xd6,0x31,0x70,0xe8,0x28,0x43,
        ], [
            0x69,0x49,0x4c,0xc6,0xad,0x47,0xe3,0xd5,0x16,0x89,0x33,0x15,0x82,0xe3,0x89,0x4e,
            0xd8,0x7a,0x73,0xc4,0xab,0x55,0x95,0x1f,0xfc,0x90,0x28,0x56,0x38,0x2d,0xbb,0x0b,
        ]),
    ];

    #[test]
    fn scale_kat() {
        let g = Point::GENERATOR;
        for (k, gx, gy) in KAT {
            let (x, y) = g.scale(&Scalar::from_u64(*k)).to_affine();
            assert_eq!(x, fe(gx), "x mismatch k={}", k);
            assert_eq!(y, fe(gy), "y mismatch k={}", k);
        }
    }

    #[test]
    fn scalar_inverse() {
        let s = Scalar::from_u64(123456789);
        assert_eq!(&s * &s.inverse(), Scalar::one());
    }
}
