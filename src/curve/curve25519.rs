//! Curve25519 / Edwards25519 defined over the prime field of order 2^255 - 19
//!
//! This module provides, for the curve commonly known as "curve25519":
//!
//! * [`FieldElement`]: an element of the base field Fp, p = 2^255 - 19
//! * [`Scalar`]: an element of the scalar field of prime order
//!   `l = 2^252 + 27742317777372353535851937790883648493`
//! * [`MontgomeryPoint`]: an x-only point on the Montgomery curve
//!   `B*y^2 = x^3 + A*x^2 + x` (with A = 486662, B = 1), with the
//!   constant-time Montgomery ladder used by X25519
//! * [`Point`]: a point on the birationally-equivalent twisted Edwards curve
//!   `-x^2 + y^2 = 1 + d*x^2*y^2` (edwards25519), with complete addition and
//!   scalar multiplication
//! * the birational map between the two forms
//!
//! The base field uses the fiat-crypto unsaturated-solinas backend, and the
//! scalar field the word-by-word Montgomery backend, in the same fashion as the
//! curves in [`crate::curve::sec2`].

use crate::curve::edwards::{EdwardsCurve, EdwardsCurveAM1};
use crate::curve::fiat::curve25519_64::*;
use crate::curve::fiat::curve25519_scalar_64::*;
use crate::curve::field::{Field, FieldSqrt, Sign};
use crate::curve::montgomery::{MontgomeryCurve, MontgomeryCurveB1};
use crate::mp::ct::{Choice, CtEqual, CtOption, CtZero};
use crate::{fiat_field_montgomery_impl, fiat_field_solinas_impl, fiat_field_sqrt_define};
#[cfg(feature = "table")]
use crate::params::curve25519::{COMB_TABLE, COMB_WINDOWS};
#[cfg(feature = "table")]
use std::convert::TryFrom;
use std::ops::{Add, Mul, Neg, Sub};

const FE_LIMBS_SIZE: usize = 5;
const GM_LIMBS_SIZE: usize = 4;

/// p = 2^255 - 19 (big-endian)
const P_BYTES: [u8; 32] = [
    0x7f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xed,
];

/// l = 2^252 + 27742317777372353535851937790883648493 (big-endian 64-bit limbs)
const ORDER_LIMBS: [u64; 4] = [
    0x1000000000000000,
    0x0000000000000000,
    0x14def9dea2f79cd6,
    0x5812631a5cf5d3ed,
];

// ***********************************************************************
// Base field Fp (unsaturated solinas) tight-domain wrappers.
//
// The fiat-crypto unsaturated-solinas routines use "loose" and "tight" field
// element bounds; the generic field macro expects operations that take and
// return tight elements, so we wrap the raw routines (relaxing inputs and
// carrying outputs) the same way `sec2::p521r1` does.
// ***********************************************************************

fn fiat_25519_nonzero(out: &mut u64, fe: &[u64; FE_LIMBS_SIZE]) {
    let mut bytes = [0u8; 32];
    let fe_tight = fiat_25519_tight_field_element(*fe);
    fiat_25519_to_bytes(&mut bytes, &fe_tight);
    *out = bytes.ct_nonzero().0;
}

const fn fiat_25519_carry_add(
    out: &mut fiat_25519_tight_field_element,
    a: &fiat_25519_tight_field_element,
    b: &fiat_25519_tight_field_element,
) {
    let mut loose = fiat_25519_loose_field_element([0u64; FE_LIMBS_SIZE]);
    fiat_25519_add(&mut loose, a, b);
    fiat_25519_carry(out, &loose);
}

const fn fiat_25519_carry_sub(
    out: &mut fiat_25519_tight_field_element,
    a: &fiat_25519_tight_field_element,
    b: &fiat_25519_tight_field_element,
) {
    let mut loose = fiat_25519_loose_field_element([0u64; FE_LIMBS_SIZE]);
    fiat_25519_sub(&mut loose, a, b);
    fiat_25519_carry(out, &loose);
}

const fn fiat_25519_mul_tight(
    out: &mut fiat_25519_tight_field_element,
    a: &fiat_25519_tight_field_element,
    b: &fiat_25519_tight_field_element,
) {
    let mut a_relaxed = fiat_25519_loose_field_element([0u64; FE_LIMBS_SIZE]);
    let mut b_relaxed = fiat_25519_loose_field_element([0u64; FE_LIMBS_SIZE]);
    fiat_25519_relax(&mut a_relaxed, a);
    fiat_25519_relax(&mut b_relaxed, b);
    fiat_25519_carry_mul(out, &a_relaxed, &b_relaxed);
}

const fn fiat_25519_square_tight(
    out: &mut fiat_25519_tight_field_element,
    a: &fiat_25519_tight_field_element,
) {
    let mut loose = fiat_25519_loose_field_element([0u64; FE_LIMBS_SIZE]);
    fiat_25519_relax(&mut loose, a);
    fiat_25519_carry_square(out, &loose);
}

const fn fiat_25519_carry_opp(
    out: &mut fiat_25519_tight_field_element,
    a: &fiat_25519_tight_field_element,
) {
    let mut loose = fiat_25519_loose_field_element([0u64; FE_LIMBS_SIZE]);
    fiat_25519_opp(&mut loose, a);
    fiat_25519_carry(out, &loose);
}

fiat_field_solinas_impl!(
    #[doc = "Element of the prime field Fp where p = 2^255 - 19"]
    FieldElement,
    255,
    P_BYTES,
    FE_LIMBS_SIZE,
    fiat_25519_tight_field_element,
    fiat_25519_loose_field_element,
    fiat_25519_carry,
    fiat_25519_relax,
    fiat_25519_nonzero,
    fiat_25519_carry_add,
    fiat_25519_carry_sub,
    fiat_25519_mul_tight,
    fiat_25519_square_tight,
    fiat_25519_carry_opp,
    fiat_25519_to_bytes,
    fiat_25519_from_bytes,
    fiat_25519_selectznz,
    le
);

/// sqrt(-1) mod p, the principal square root of -1 (big-endian)
const SQRT_M1_BYTES: [u8; 32] = [
    0x2b, 0x83, 0x24, 0x80, 0x4f, 0xc1, 0xdf, 0x0b, 0x2b, 0x4d, 0x00, 0x99, 0x3d, 0xfb, 0xd7, 0xa7,
    0x2f, 0x43, 0x18, 0x06, 0xad, 0x2f, 0xe4, 0x78, 0xc4, 0xee, 0x1b, 0x27, 0x4a, 0x0e, 0xa0, 0xb0,
];

impl FieldElement {
    const SQRT_M1: FieldElement = FieldElement::from_bytes_unchecked_be(&SQRT_M1_BYTES);

    /// Returns `(self^(2^250 - 1), self^11)`.
    ///
    /// This is the shared prefix of both the inverse (`self^(p-2)`) and the
    /// `(p-5)/8` exponentiation used by the square root, expressed with the
    /// classic ref10 addition chain.
    fn pow_2_250_m1(&self) -> (FieldElement, FieldElement) {
        let z1 = self;
        let z2 = z1.square();
        let z4 = z2.square();
        let z8 = z4.square();
        let z9 = &z8 * z1;
        let z11 = &z9 * &z2;
        let z22 = z11.square();
        let z_5_0 = &z22 * &z9;
        let t = z_5_0.square_rep(5);
        let z_10_0 = &t * &z_5_0;
        let t = z_10_0.square_rep(10);
        let z_20_0 = &t * &z_10_0;
        let t = z_20_0.square_rep(20);
        let t = &t * &z_20_0;
        let t = t.square_rep(10);
        let z_50_0 = &t * &z_10_0;
        let t = z_50_0.square_rep(50);
        let z_100_0 = &t * &z_50_0;
        let t = z_100_0.square_rep(100);
        let t = &t * &z_100_0;
        let t = t.square_rep(50);
        let t = &t * &z_50_0;
        (t, z11)
    }

    /// Multiplicative inverse without the non-zero precondition: returns 0 for 0.
    fn invert_or_zero(&self) -> Self {
        // p - 2 = 2^255 - 21
        let (t, z11) = self.pow_2_250_m1();
        &t.square_rep(5) * &z11
    }

    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    pub fn inverse(&self) -> Self {
        assert!(!self.is_zero());
        self.invert_or_zero()
    }

    /// Compute a square root 'x' of the field element such that x*x = self
    ///
    /// As p ≡ 5 (mod 8), this uses the Atkin method: the candidate root is
    /// `self^((p+3)/8)`, corrected by `sqrt(-1)` when it is a root of `-self`
    /// instead. `None` is returned when `self` is not a quadratic residue.
    pub fn sqrt(&self) -> CtOption<Self> {
        let (t, _z11) = self.pow_2_250_m1();
        // self^(2^252 - 3) = (2^250 - 1) doubled twice, times self
        let p58 = &t.square_rep(2) * self;
        // candidate = self^((p+3)/8) = self^(2^252 - 2)
        let cand = &p58 * self;
        let cand2 = &cand * &Self::SQRT_M1;

        let is_root = cand.square().ct_eq(self);
        let r = Self::ct_select(is_root, &cand, &cand2);
        let present = r.square().ct_eq(self);
        CtOption::from((present, r))
    }
}
fiat_field_sqrt_define!(FieldElement);

fiat_field_montgomery_impl!(
    #[doc = "Element of the scalar field of prime order l = 2^252 + 27742317777372353535851937790883648493"]
    Scalar,
    253,
    ORDER_LIMBS,
    GM_LIMBS_SIZE,
    fiat_25519_scalar_non_montgomery_domain_field_element,
    fiat_25519_scalar_nonzero,
    fiat_25519_scalar_add,
    fiat_25519_scalar_sub,
    fiat_25519_scalar_mul,
    fiat_25519_scalar_square,
    fiat_25519_scalar_opp,
    fiat_25519_scalar_to_bytes,
    fiat_25519_scalar_from_bytes,
    fiat_25519_scalar_montgomery_domain_field_element,
    fiat_25519_scalar_to_montgomery,
    fiat_25519_scalar_from_montgomery,
    fiat_25519_scalar_selectznz,
    fiat_25519_scalar_msat,
    fiat_25519_scalar_divstep,
    fiat_25519_scalar_divstep_precomp,
    le
);

// l - 2, big-endian, for the Fermat scalar inversion
const ORDER_M2_BYTES: [u8; 32] = [
    0x10, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x14, 0xde, 0xf9, 0xde, 0xa2, 0xf7, 0x9c, 0xd6, 0x58, 0x12, 0x63, 0x1a, 0x5c, 0xf5, 0xd3, 0xeb,
];

impl Scalar {
    /// Get the multiplicative inverse, computed as `self^(l-2)` (Fermat) using a
    /// fixed 4-bit window over the constant exponent `l-2`.
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic.
    pub fn inverse(&self) -> Self {
        assert!(!self.is_zero());

        // table[i] = self^i for i in 0..16
        let mut table: [Scalar; 16] = core::array::from_fn(|_| Scalar::one());
        table[1] = self.clone();
        for i in 2..16 {
            table[i] = &table[i - 1] * self;
        }

        // MSB-first windowed square-and-multiply over the bytes of l-2; the
        // exponent is a public constant so the control flow leaks nothing.
        let mut acc = Scalar::one();
        let mut started = false;
        for &byte in ORDER_M2_BYTES.iter() {
            for &nibble in &[byte >> 4, byte & 0x0f] {
                if started {
                    acc = acc.square_rep(4);
                }
                if nibble != 0 {
                    acc = if started {
                        &acc * &table[nibble as usize]
                    } else {
                        table[nibble as usize].clone()
                    };
                    started = true;
                }
            }
        }
        acc
    }
}

// ***********************************************************************
// Curve definitions
// ***********************************************************************

// Montgomery form: B y^2 = x^3 + A x^2 + x , with A = 486662, B = 1
const MONT_A_BYTES: [u8; 32] = [
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x07, 0x6d, 0x06,
];
const MONT_B_BYTES: [u8; 32] = [
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
];
// (A + 2) / 4 = 121666
const MONT_A24_BYTES: [u8; 32] = [
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0xdb, 0x42,
];
// Montgomery base point u-coordinate = 9
const MONT_GU_BYTES: [u8; 32] = [
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x09,
];
// Montgomery base point v-coordinate (used for the birational map cross-check)
#[allow(dead_code)]
const MONT_GV_BYTES: [u8; 32] = [
    0x20, 0xae, 0x19, 0xa1, 0xb8, 0xa0, 0x86, 0xb4, 0xe0, 0x1e, 0xdd, 0x2c, 0x77, 0x48, 0xd1, 0x4c,
    0x92, 0x3d, 0x4d, 0x7e, 0x6d, 0x7c, 0x61, 0xb2, 0x29, 0xe9, 0xc5, 0xa2, 0x7e, 0xce, 0xd3, 0xd9,
];

// Twisted Edwards form: a x^2 + y^2 = 1 + d x^2 y^2 , with a = -1
const ED_A_BYTES: [u8; 32] = [
    0x7f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xec,
];
const ED_D_BYTES: [u8; 32] = [
    0x52, 0x03, 0x6c, 0xee, 0x2b, 0x6f, 0xfe, 0x73, 0x8c, 0xc7, 0x40, 0x79, 0x77, 0x79, 0xe8, 0x98,
    0x00, 0x70, 0x0a, 0x4d, 0x41, 0x41, 0xd8, 0xab, 0x75, 0xeb, 0x4d, 0xca, 0x13, 0x59, 0x78, 0xa3,
];
// 2 * d
const ED_D2_BYTES: [u8; 32] = [
    0x24, 0x06, 0xd9, 0xdc, 0x56, 0xdf, 0xfc, 0xe7, 0x19, 0x8e, 0x80, 0xf2, 0xee, 0xf3, 0xd1, 0x30,
    0x00, 0xe0, 0x14, 0x9a, 0x82, 0x83, 0xb1, 0x56, 0xeb, 0xd6, 0x9b, 0x94, 0x26, 0xb2, 0xf1, 0x59,
];
// Edwards base point coordinates
const ED_GX_BYTES: [u8; 32] = [
    0x21, 0x69, 0x36, 0xd3, 0xcd, 0x6e, 0x53, 0xfe, 0xc0, 0xa4, 0xe2, 0x31, 0xfd, 0xd6, 0xdc, 0x5c,
    0x69, 0x2c, 0xc7, 0x60, 0x95, 0x25, 0xa7, 0xb2, 0xc9, 0x56, 0x2d, 0x60, 0x8f, 0x25, 0xd5, 0x1a,
];
const ED_GY_BYTES: [u8; 32] = [
    0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
    0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x58,
];
// Edwards base point extended T = GX * GY
const ED_GT_BYTES: [u8; 32] = [
    0x67, 0x87, 0x5f, 0x0f, 0xd7, 0x8b, 0x76, 0x65, 0x66, 0xea, 0x4e, 0x8e, 0x64, 0xab, 0xe3, 0x7d,
    0x20, 0xf0, 0x9f, 0x80, 0x77, 0x51, 0x52, 0xf5, 0x6d, 0xde, 0x8a, 0xb3, 0xa5, 0xb7, 0xdd, 0xa3,
];

// c = sqrt(-486664), the constant relating the Montgomery and Edwards forms
const MAP_C_BYTES: [u8; 32] = [
    0x70, 0xd9, 0x12, 0x0b, 0x9f, 0x5f, 0xf9, 0x44, 0x2d, 0x84, 0xf7, 0x23, 0xfc, 0x03, 0xb0, 0x81,
    0x3a, 0x5e, 0x2c, 0x2e, 0xb4, 0x82, 0xe5, 0x7d, 0x33, 0x91, 0xfb, 0x55, 0x00, 0xba, 0x81, 0xe7,
];
const MAP_C: FieldElement = FieldElement::from_bytes_unchecked_be(&MAP_C_BYTES);

/// The Montgomery curve `y^2 = x^3 + 486662 x^2 + x` (curve25519)
#[derive(Debug, Clone, Copy)]
pub struct Curve;

impl MontgomeryCurve for Curve {
    type FieldElement = FieldElement;
    const A: FieldElement = FieldElement::from_bytes_unchecked_be(&MONT_A_BYTES);
    const B: FieldElement = FieldElement::from_bytes_unchecked_be(&MONT_B_BYTES);
    const A24: FieldElement = FieldElement::from_bytes_unchecked_be(&MONT_A24_BYTES);
}
impl MontgomeryCurveB1 for Curve {}

/// The twisted Edwards curve `-x^2 + y^2 = 1 + d*x^2*y^2` (edwards25519)
#[derive(Debug, Clone, Copy)]
pub struct EdCurve;

impl EdwardsCurve for EdCurve {
    type FieldElement = FieldElement;
    const A: FieldElement = FieldElement::from_bytes_unchecked_be(&ED_A_BYTES);
    const D: FieldElement = FieldElement::from_bytes_unchecked_be(&ED_D_BYTES);
    const D2: FieldElement = FieldElement::from_bytes_unchecked_be(&ED_D2_BYTES);
}
impl EdwardsCurveAM1 for EdCurve {}

// ***********************************************************************
// Montgomery x-only point and ladder (X25519)
// ***********************************************************************

/// x-only point on the Montgomery curve, in projective `(X:Z)` coordinates.
///
/// Only the x-coordinate is tracked; this representation supports the
/// differential (Montgomery) ladder but not general point addition.
#[derive(Clone, Debug)]
pub struct MontgomeryPoint {
    x: FieldElement,
    z: FieldElement,
}

/// Constant-time conditional swap of two field elements.
fn cswap(swap: Choice, a: &mut FieldElement, b: &mut FieldElement) {
    let na = FieldElement::ct_select(swap, b, a);
    let nb = FieldElement::ct_select(swap, a, b);
    *a = na;
    *b = nb;
}

/// Constant-time Montgomery ladder: returns the u-coordinate of `k * P`, where
/// `P` is the affine u-coordinate `base_u` and `k` is the big-endian scalar.
fn ladder(base_u: &FieldElement, k_be: &[u8]) -> FieldElement {
    let a24 = <Curve as MontgomeryCurve>::A24;
    let x1 = base_u.clone();
    let mut x2 = FieldElement::one();
    let mut z2 = FieldElement::zero();
    let mut x3 = base_u.clone();
    let mut z3 = FieldElement::one();
    let mut swap = 0u8;

    for &byte in k_be.iter() {
        for i in (0..8).rev() {
            let bit = (byte >> i) & 1;
            swap ^= bit;
            let sw = Choice((swap as u64) & 1);
            cswap(sw, &mut x2, &mut x3);
            cswap(sw, &mut z2, &mut z3);
            swap = bit;

            // differential add-and-double
            let a = &x2 + &z2;
            let aa = a.square();
            let b = &x2 - &z2;
            let bb = b.square();
            let e = &aa - &bb;
            let c = &x3 + &z3;
            let d = &x3 - &z3;
            let da = &d * &a;
            let cb = &c * &b;
            x3 = (&da + &cb).square();
            z3 = &x1 * &(&da - &cb).square();
            x2 = &aa * &bb;
            z2 = &e * &(&bb + &(&a24 * &e));
        }
    }
    let sw = Choice((swap as u64) & 1);
    cswap(sw, &mut x2, &mut x3);
    cswap(sw, &mut z2, &mut z3);

    &x2 * &z2.invert_or_zero()
}

impl MontgomeryPoint {
    /// Base point with u-coordinate 9
    pub const GENERATOR: Self = MontgomeryPoint {
        x: FieldElement::from_bytes_unchecked_be(&MONT_GU_BYTES),
        z: FieldElement::one(),
    };

    /// Build an x-only point from its affine u-coordinate
    pub fn from_u(u: &FieldElement) -> Self {
        MontgomeryPoint {
            x: u.clone(),
            z: FieldElement::one(),
        }
    }

    /// Return the affine u-coordinate (0 if this is the point at infinity)
    pub fn u(&self) -> FieldElement {
        &self.x * &self.z.invert_or_zero()
    }

    /// Scalar multiplication by a big-endian scalar using the Montgomery ladder.
    pub fn scale_bytes(&self, k_be: &[u8]) -> Self {
        MontgomeryPoint {
            x: ladder(&self.u(), k_be),
            z: FieldElement::one(),
        }
    }

    /// Scalar multiplication by a [`Scalar`] using the Montgomery ladder.
    pub fn scale(&self, k: &Scalar) -> Self {
        self.scale_bytes(&k.to_bytes_be())
    }

    /// Map this Montgomery point to a twisted Edwards point, picking the
    /// `x_sign` branch of the recovered x-coordinate.
    ///
    /// Returns `None` for the points where the rational map is undefined.
    pub fn to_edwards(&self, x_sign: Sign) -> Option<Point> {
        let u = self.u();
        let one = FieldElement::one();
        let denom = &u + &one;
        if denom.is_zero() {
            return None;
        }
        // y = (u - 1) / (u + 1)
        let y = &(&u - &one) * &denom.inverse();
        Point::decompress(&y, x_sign)
    }
}

impl PartialEq for MontgomeryPoint {
    fn eq(&self, other: &Self) -> bool {
        // (X1:Z1) == (X2:Z2) iff X1*Z2 == X2*Z1
        let l = &self.x * &other.z;
        let r = &other.x * &self.z;
        l.ct_eq(&r).is_true()
    }
}
impl Eq for MontgomeryPoint {}

impl<'a, 'b> Mul<&'b Scalar> for &'a MontgomeryPoint {
    type Output = MontgomeryPoint;
    fn mul(self, k: &'b Scalar) -> MontgomeryPoint {
        self.scale(k)
    }
}

// ***********************************************************************
// Twisted Edwards point (full group)
// ***********************************************************************

/// Point on the twisted Edwards curve edwards25519, in extended homogeneous
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

    /// Curve generator (base point)
    pub const GENERATOR: Self = Point {
        x: FieldElement::from_bytes_unchecked_be(&ED_GX_BYTES),
        y: FieldElement::from_bytes_unchecked_be(&ED_GY_BYTES),
        z: FieldElement::one(),
        t: FieldElement::from_bytes_unchecked_be(&ED_GT_BYTES),
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

    fn ct_select(cond: Choice, a: &Point, b: &Point) -> Point {
        Point {
            x: FieldElement::ct_select(cond, &a.x, &b.x),
            y: FieldElement::ct_select(cond, &a.y, &b.y),
            z: FieldElement::ct_select(cond, &a.z, &b.z),
            t: FieldElement::ct_select(cond, &a.t, &b.t),
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
                let added = q.add(self);
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

    /// Map this Edwards point to the Montgomery u-coordinate `(1 + y)/(1 - y)`.
    ///
    /// Returns `None` when `y == 1` (the rational map sends it to infinity).
    pub fn to_montgomery_u(&self) -> Option<FieldElement> {
        let (_x, y) = self.to_affine();
        let one = FieldElement::one();
        let denom = &one - &y;
        if denom.is_zero() {
            return None;
        }
        Some(&(&one + &y) * &denom.inverse())
    }

    /// Map this Edwards point to the full Montgomery coordinates `(u, v)`.
    ///
    /// Returns `None` for the points where the rational map is undefined.
    pub fn to_montgomery(&self) -> Option<(FieldElement, FieldElement)> {
        let (x, y) = self.to_affine();
        if x.is_zero() {
            return None;
        }
        let one = FieldElement::one();
        let denom = &one - &y;
        if denom.is_zero() {
            return None;
        }
        let u = &(&one + &y) * &denom.inverse();
        // v = c*u/x
        let v = &(&MAP_C * &u) * &x.invert_or_zero();
        Some((u, v))
    }

    /// Build an Edwards point from full Montgomery coordinates `(u, v)`.
    pub fn from_montgomery(u: &FieldElement, v: &FieldElement) -> Option<Self> {
        let one = FieldElement::one();
        let denom = u + &one;
        if denom.is_zero() || v.is_zero() {
            return None;
        }
        // y = (u - 1)/(u + 1) ; x = c*u/v
        let y = &(u - &one) * &denom.inverse();
        let x = &(&MAP_C * u) * &v.inverse();
        Point::from_coordinate(&x, &y)
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
    static V: std::sync::OnceLock<Box<[[Point; 16]; COMB_WINDOWS]>> = std::sync::OnceLock::new();
    &**V.get_or_init(build_comb_table)
}

#[cfg(feature = "table")]
fn build_comb_table() -> Box<[[Point; 16]; COMB_WINDOWS]> {
    let mut windows: Vec<[Point; 16]> = Vec::with_capacity(COMB_WINDOWS);
    for row in COMB_TABLE.iter() {
        let mut window: [Point; 16] = core::array::from_fn(|_| Point::IDENTITY);
        for (slot, (x, y)) in window.iter_mut().skip(1).zip(row.iter()) {
            *slot = Point::from_affine(
                &FieldElement::from_bytes_unchecked_le(x),
                &FieldElement::from_bytes_unchecked_le(y),
            );
        }
        windows.push(window);
    }
    <Box<[[Point; 16]; COMB_WINDOWS]>>::try_from(windows.into_boxed_slice())
        .ok()
        .expect("comb window count matches COMB_WINDOWS")
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
impl<'a, 'b> Sub<&'b Point> for &'a Point {
    type Output = Point;
    fn sub(self, other: &'b Point) -> Point {
        Point::add(self, &other.negate())
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a Point {
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

#[cfg(test)]
mod tests {
    mod fe {
        use super::super::FieldElement;
        use crate::{fiat_field_sqrt_unittest, fiat_field_unittest};
        fiat_field_unittest!(FieldElement);
        fiat_field_sqrt_unittest!(FieldElement);
    }
    mod sc {
        use super::super::Scalar;
        use crate::fiat_field_unittest;
        fiat_field_unittest!(Scalar);
        crate::fiat_field_safegcd_unittest!(Scalar);
    }
    mod curve {
        use super::super::*;

        // l, the prime order of the base point (big-endian)
        const ORDER_BYTES: [u8; 32] = [
            0x10, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x14, 0xde, 0xf9, 0xde, 0xa2, 0xf7, 0x9c, 0xd6, 0x58, 0x12, 0x63, 0x1a,
            0x5c, 0xf5, 0xd3, 0xed,
        ];

        // the in-module test vectors are big-endian
        fn fe(bytes: &[u8; 32]) -> FieldElement {
            FieldElement::from_bytes_be(bytes).expect("valid field element")
        }

        #[test]
        fn edwards_generator_on_curve() {
            let (x, y) = Point::GENERATOR.to_affine();
            assert!(Point::from_coordinate(&x, &y).is_some());
            // the stored extended T must equal x*y
            assert_eq!(&x * &y, fe(&super::super::ED_GT_BYTES));
        }

        #[test]
        fn edwards_identity() {
            let g = Point::GENERATOR;
            assert_eq!(&g + &Point::IDENTITY, g);
            assert_eq!(&Point::IDENTITY + &g, g);
            assert_eq!(&g + &(-&g), Point::IDENTITY);
        }

        #[test]
        fn edwards_double_matches_add() {
            let g = Point::GENERATOR;
            assert_eq!(g.double(), &g + &g);
            assert_eq!(g.scale_bytes(&[2]), g.double());
            assert_eq!(g.scale_bytes(&[3]), &g.double() + &g);
        }

        #[test]
        fn edwards_generator_has_order_l() {
            // l * G == identity, and (l - 1)*G == -G
            assert_eq!(Point::GENERATOR.scale_bytes(&ORDER_BYTES), Point::IDENTITY);
        }

        #[test]
        fn mul_base_matches_scale() {
            // the fixed-base comb must agree with the generic scalar mult
            assert_eq!(Point::mul_base(&Scalar::from_u64(0)), Point::IDENTITY);
            for k in [1u64, 2, 9, 16, 255, 256, 1234567, 0x0123_4567_89ab_cdef] {
                let s = Scalar::from_u64(k);
                assert_eq!(Point::mul_base(&s), Point::GENERATOR.scale(&s), "k={}", k);
            }
            // a full-width scalar too
            let mut s = Scalar::from_u64(0x0123_4567_89ab_cdef);
            for _ in 0..5 {
                s = s.square();
            }
            assert_eq!(Point::mul_base(&s), Point::GENERATOR.scale(&s));
        }

        #[test]
        fn edwards_scalar_linearity() {
            let g = Point::GENERATOR;
            for k in [1u8, 2, 5, 9, 200] {
                let lhs = g.scale_bytes(&[k]);
                // add G to itself k times
                let mut rhs = Point::IDENTITY;
                for _ in 0..k {
                    rhs = &rhs + &g;
                }
                assert_eq!(lhs, rhs, "k={}", k);
            }
        }

        #[test]
        fn edwards_compress_roundtrip() {
            let g = Point::GENERATOR;
            for k in [1u8, 2, 3, 7, 50] {
                let p = g.scale_bytes(&[k]);
                let (y, sign) = p.compress();
                let q = Point::decompress(&y, sign).expect("decompress");
                assert_eq!(p, q, "k={}", k);
            }
        }

        // independently (python) computed u-coordinates of k*(u=9) on the
        // Montgomery curve
        const LADDER_9_2: [u8; 32] = [
            0x20, 0xd3, 0x42, 0xd5, 0x18, 0x73, 0xf1, 0xb7, 0xd9, 0x75, 0x0c, 0x68, 0x7d, 0x15,
            0x71, 0x14, 0x8f, 0x3f, 0x5c, 0xed, 0x1e, 0x35, 0x0b, 0x5c, 0x5c, 0xae, 0x46, 0x9c,
            0xdd, 0x68, 0x4e, 0xfb,
        ];
        const LADDER_9_5: [u8; 32] = [
            0x41, 0xb6, 0xec, 0x3c, 0x50, 0xee, 0x7a, 0xf2, 0x03, 0xc0, 0x02, 0x6e, 0x5e, 0x07,
            0x9e, 0x7f, 0xa8, 0xcb, 0xc9, 0xbc, 0x58, 0x1d, 0x49, 0xcb, 0x0d, 0x53, 0x7d, 0x57,
            0x78, 0x49, 0x7c, 0x87,
        ];
        const LADDER_9_7: [u8; 32] = [
            0x0d, 0xaf, 0x32, 0xe7, 0xed, 0x80, 0x99, 0x12, 0x2b, 0x2d, 0xfa, 0x4c, 0x1d, 0x8c,
            0x4a, 0x20, 0xc0, 0x97, 0x2a, 0x15, 0x38, 0xbf, 0x05, 0x75, 0x33, 0x8a, 0xae, 0x0f,
            0xe0, 0x84, 0x18, 0x28,
        ];

        #[test]
        fn montgomery_ladder_kat() {
            let g = MontgomeryPoint::GENERATOR;
            assert_eq!(g.scale_bytes(&[1]).u(), fe(&super::super::MONT_GU_BYTES));
            assert_eq!(g.scale_bytes(&[2]).u(), fe(&LADDER_9_2));
            assert_eq!(g.scale_bytes(&[5]).u(), fe(&LADDER_9_5));
            assert_eq!(g.scale_bytes(&[7]).u(), fe(&LADDER_9_7));
        }

        #[test]
        fn ladder_matches_edwards() {
            // [k](u=9) via the ladder == montgomery-u of [k]*G_edwards
            let mg = MontgomeryPoint::GENERATOR;
            let eg = Point::GENERATOR;
            for k in [1u64, 2, 3, 4, 5, 9, 1000, 1234567] {
                let s = Scalar::from_u64(k);
                let via_ladder = mg.scale(&s).u();
                let via_edwards = eg.scale(&s).to_montgomery_u().expect("u defined");
                assert_eq!(via_ladder, via_edwards, "k={}", k);
            }
        }

        #[test]
        fn map_generator() {
            // Edwards generator maps to the Montgomery base point (9, v)
            let (u, v) = Point::GENERATOR.to_montgomery().expect("map defined");
            assert_eq!(u, fe(&super::super::MONT_GU_BYTES));
            assert_eq!(v, fe(&super::super::MONT_GV_BYTES));
            // and the inverse map recovers the generator
            let back = Point::from_montgomery(&u, &v).expect("inverse map");
            assert_eq!(back, Point::GENERATOR);
        }

        #[test]
        fn scalar_inverse() {
            let s = Scalar::from_u64(123456789);
            assert_eq!(&s * &s.inverse(), Scalar::one());
        }

        #[test]
        fn default_endianness_is_le() {
            // curve25519's field and scalar default to little-endian
            let x = fe(&super::super::ED_GX_BYTES);
            assert_eq!(x.to_bytes(), x.to_bytes_le());
            assert_ne!(x.to_bytes(), x.to_bytes_be());
            assert_eq!(FieldElement::from_bytes(&x.to_bytes_le()).unwrap(), x);

            let s = Scalar::from_u64(0x1234_5678_9abc_def0);
            assert_eq!(s.to_bytes(), s.to_bytes_le());
            assert_eq!(Scalar::from_bytes(&s.to_bytes_le()).unwrap(), s);
        }
    }
}
