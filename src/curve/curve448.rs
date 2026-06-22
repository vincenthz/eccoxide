//! Curve448 ("Goldilocks") defined over the prime field of order 2^448 - 2^224 - 1
//!
//! This module provides, for the curve commonly known as "curve448":
//!
//! * [`FieldElement`]: an element of the base field Fp, p = 2^448 - 2^224 - 1
//! * [`MontgomeryPoint`]: an x-only point on the Montgomery curve
//!   `B·y^2 = x^3 + A·x^2 + x` (with A = 156326, B = 1), with the
//!   constant-time Montgomery ladder used by X448
//!
//! The base field uses the fiat-crypto unsaturated-solinas backend, in the same
//! fashion as the curves in [`crate::curve::sec2`]. Like curve25519, the wire
//! byte order is little-endian.
//!
//! Note: there is no scalar field type here. X448 multiplies by the clamped
//! scalar integer directly (not reduced modulo the group order), so it consumes
//! raw scalar bytes; a reduced [`Scalar`](crate::curve::curve25519::Scalar)-style
//! type would only be needed for a signature scheme (Ed448), which is not
//! provided.

use crate::curve::field::{Field, FieldSqrt, Sign};
use crate::curve::fiat::p448_solinas_64::*;
use crate::curve::montgomery::{MontgomeryCurve, MontgomeryCurveB1};
use crate::mp::ct::{Choice, CtEqual, CtOption, CtSelect, CtZero};
use crate::{fiat_field_solinas_impl, fiat_field_sqrt_define};

const FE_LIMBS_SIZE: usize = 8;

/// p = 2^448 - 2^224 - 1 (big-endian)
const P_BYTES: [u8; 56] = [
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
];

// ***********************************************************************
// Base field Fp (unsaturated solinas) tight-domain wrappers.
//
// The fiat-crypto unsaturated-solinas routines use "loose" and "tight" field
// element bounds; the generic field macro expects operations that take and
// return tight elements, so we wrap the raw routines (relaxing inputs and
// carrying outputs) the same way `sec2::p521r1` and `curve25519` do.
// ***********************************************************************

fn fiat_p448_nonzero(out: &mut u64, fe: &[u64; FE_LIMBS_SIZE]) {
    let mut bytes = [0u8; 56];
    let fe_tight = fiat_p448_tight_field_element(*fe);
    fiat_p448_to_bytes(&mut bytes, &fe_tight);
    *out = bytes.ct_nonzero().0;
}

const fn fiat_p448_carry_add(
    out: &mut fiat_p448_tight_field_element,
    a: &fiat_p448_tight_field_element,
    b: &fiat_p448_tight_field_element,
) {
    let mut loose = fiat_p448_loose_field_element([0u64; FE_LIMBS_SIZE]);
    fiat_p448_add(&mut loose, a, b);
    fiat_p448_carry(out, &loose);
}

const fn fiat_p448_carry_sub(
    out: &mut fiat_p448_tight_field_element,
    a: &fiat_p448_tight_field_element,
    b: &fiat_p448_tight_field_element,
) {
    let mut loose = fiat_p448_loose_field_element([0u64; FE_LIMBS_SIZE]);
    fiat_p448_sub(&mut loose, a, b);
    fiat_p448_carry(out, &loose);
}

const fn fiat_p448_mul_tight(
    out: &mut fiat_p448_tight_field_element,
    a: &fiat_p448_tight_field_element,
    b: &fiat_p448_tight_field_element,
) {
    let mut a_relaxed = fiat_p448_loose_field_element([0u64; FE_LIMBS_SIZE]);
    let mut b_relaxed = fiat_p448_loose_field_element([0u64; FE_LIMBS_SIZE]);
    fiat_p448_relax(&mut a_relaxed, a);
    fiat_p448_relax(&mut b_relaxed, b);
    fiat_p448_carry_mul(out, &a_relaxed, &b_relaxed);
}

const fn fiat_p448_square_tight(
    out: &mut fiat_p448_tight_field_element,
    a: &fiat_p448_tight_field_element,
) {
    let mut loose = fiat_p448_loose_field_element([0u64; FE_LIMBS_SIZE]);
    fiat_p448_relax(&mut loose, a);
    fiat_p448_carry_square(out, &loose);
}

const fn fiat_p448_carry_opp(
    out: &mut fiat_p448_tight_field_element,
    a: &fiat_p448_tight_field_element,
) {
    let mut loose = fiat_p448_loose_field_element([0u64; FE_LIMBS_SIZE]);
    fiat_p448_opp(&mut loose, a);
    fiat_p448_carry(out, &loose);
}

fiat_field_solinas_impl!(
    #[doc = "Element of the prime field Fp where p = 2^448 - 2^224 - 1"]
    FieldElement,
    448,
    P_BYTES,
    FE_LIMBS_SIZE,
    fiat_p448_tight_field_element,
    fiat_p448_loose_field_element,
    fiat_p448_carry,
    fiat_p448_relax,
    fiat_p448_nonzero,
    fiat_p448_carry_add,
    fiat_p448_carry_sub,
    fiat_p448_mul_tight,
    fiat_p448_square_tight,
    fiat_p448_carry_opp,
    fiat_p448_to_bytes,
    fiat_p448_from_bytes,
    fiat_p448_selectznz,
    le
);

/// p - 2 (big-endian), the exponent for the Fermat inverse a^(p-2)
const PM2_BYTES: [u8; 56] = [
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfd,
];

/// (p + 1) / 4 (big-endian), the exponent for the square root (p ≡ 3 mod 4)
const P14_BYTES: [u8; 56] = [
    0x3f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xc0, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

impl FieldElement {
    /// Fixed 4-bit windowed exponentiation by the big-endian constant exponent
    /// `exp`. The exponent is public, so the control flow leaks nothing.
    fn pow_const(&self, exp: &[u8]) -> Self {
        // table[i] = self^i for i in 0..16
        let mut table: [FieldElement; 16] = core::array::from_fn(|_| FieldElement::one());
        table[1] = self.clone();
        for i in 2..16 {
            table[i] = &table[i - 1] * self;
        }

        let mut acc = FieldElement::one();
        let mut started = false;
        for &byte in exp.iter() {
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

    /// Multiplicative inverse without the non-zero precondition: returns 0 for 0.
    fn invert_or_zero(&self) -> Self {
        self.pow_const(&PM2_BYTES)
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
    /// As p ≡ 3 (mod 4), the candidate root is `self^((p+1)/4)`. `None` is
    /// returned when `self` is not a quadratic residue.
    pub fn sqrt(&self) -> CtOption<Self> {
        let r = self.pow_const(&P14_BYTES);
        let r2 = &r * &r;
        CtOption::from((r2.ct_eq(self), r))
    }
}
fiat_field_sqrt_define!(FieldElement);

// ***********************************************************************
// Curve448 Montgomery form: B y^2 = x^3 + A x^2 + x , with A = 156326, B = 1
// ***********************************************************************

const MONT_A_BYTES: [u8; 56] = [
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x62, 0xa6,
];
const MONT_B_BYTES: [u8; 56] = [
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
];
// (A + 2) / 4 = 39082
const MONT_A24_BYTES: [u8; 56] = [
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x98, 0xaa,
];
// Curve448 base point u-coordinate = 5
const MONT_GU_BYTES: [u8; 56] = [
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x05,
];

/// The Montgomery curve `y^2 = x^3 + 156326·x^2 + x` (curve448)
#[derive(Debug, Clone, Copy)]
pub struct Curve;

impl MontgomeryCurve for Curve {
    type FieldElement = FieldElement;
    const A: FieldElement = FieldElement::from_bytes_unchecked_be(&MONT_A_BYTES);
    const B: FieldElement = FieldElement::from_bytes_unchecked_be(&MONT_B_BYTES);
    const A24: FieldElement = FieldElement::from_bytes_unchecked_be(&MONT_A24_BYTES);
}
impl MontgomeryCurveB1 for Curve {}

// ***********************************************************************
// Montgomery x-only point and ladder (X448)
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

/// Constant-time Montgomery ladder: returns the u-coordinate of `k · P`, where
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
    /// Base point with u-coordinate 5
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
}

impl PartialEq for MontgomeryPoint {
    fn eq(&self, other: &Self) -> bool {
        // (X1:Z1) == (X2:Z2) iff X1·Z2 == X2·Z1
        let l = &self.x * &other.z;
        let r = &other.x * &self.z;
        l.ct_eq(&r).is_true()
    }
}
impl Eq for MontgomeryPoint {}

#[cfg(test)]
mod tests {
    mod fe {
        use super::super::FieldElement;
        use crate::{fiat_field_sqrt_unittest, fiat_field_unittest};
        fiat_field_unittest!(FieldElement);
        fiat_field_sqrt_unittest!(FieldElement);
    }

    mod curve {
        use super::super::*;

        fn fe_be(bytes: &[u8; 56]) -> FieldElement {
            FieldElement::from_bytes_be(bytes).expect("valid field element")
        }

        #[test]
        fn base_point_on_curve() {
            // v^2 = u^3 + A·u^2 + u must be a quadratic residue for u = 5
            let u = fe_be(&super::super::MONT_GU_BYTES);
            let a = <Curve as MontgomeryCurve>::A;
            let vv = &(&(&u.square() * &u) + &(&a * &u.square())) + &u;
            assert!(vv.sqrt().into_option().is_some(), "base u=5 not on curve");
        }

        #[test]
        fn ladder_by_one_is_identity_on_u() {
            let g = MontgomeryPoint::GENERATOR;
            assert_eq!(g.scale_bytes(&[1]).u(), fe_be(&super::super::MONT_GU_BYTES));
        }

        #[test]
        fn ladder_diffie_hellman_commutes() {
            // [a]([b]G) == [b]([a]G) — the core DH property of the ladder
            let g = MontgomeryPoint::GENERATOR;
            let a: [u8; 56] = core::array::from_fn(|i| (i as u8).wrapping_mul(7).wrapping_add(3));
            let b: [u8; 56] = core::array::from_fn(|i| (i as u8).wrapping_mul(5).wrapping_add(9));
            let ab = MontgomeryPoint::from_u(&g.scale_bytes(&b).u()).scale_bytes(&a);
            let ba = MontgomeryPoint::from_u(&g.scale_bytes(&a).u()).scale_bytes(&b);
            assert_eq!(ab.u(), ba.u());
        }

        #[test]
        fn field_inverse_and_sqrt() {
            let x = FieldElement::from_u64(123456789);
            assert_eq!(&x * &x.inverse(), FieldElement::one());
            let sq = x.square();
            let r = sq.sqrt().into_option().expect("square is a residue");
            assert_eq!(&r * &r, sq);
        }
    }
}
