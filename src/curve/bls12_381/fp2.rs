//! BLS12-381 quadratic extension field `Fp2`
//!
//! `Fp2 = Fp[u]/(u^2 + 1)`, i.e. elements are `c0 + c1*u` with `u^2 = -1`
//! (`p = 3 (mod 4)` guarantees `x^2 + 1` is irreducible over [`Fp`]). This is
//! the field the G2 twist is defined over.
//!
//! It is implemented on top of [`Fp`] and wired into the generic [`Field`] /
//! [`FieldSqrt`] abstractions so the shared curve machinery (affine /
//! projective points, complete addition, comb) works over it unchanged.

use super::fp::Fp;
use crate::curve::field::{Field, FieldSqrt, Sign};
use crate::mp::ct::{Choice, CtEqual, CtOption, CtSelect, CtZero};
use crate::params::bls12_381::{P_MINUS1_DIV2_BYTES, P_MINUS3_DIV4_BYTES};
use std::ops::{Add, Mul, Neg, Sub};

/// Element `c0 + c1*u` of the quadratic extension field `Fp2 = Fp[u]/(u^2 + 1)`.
#[derive(Clone)]
pub struct Fp2 {
    pub c0: Fp,
    pub c1: Fp,
}

impl Fp2 {
    /// Size in bytes of the serialized element (two [`Fp`] elements).
    pub const SIZE_BYTES: usize = 2 * Fp::SIZE_BYTES;

    /// Build an element from its two components `c0 + c1*u`.
    pub const fn new(c0: Fp, c1: Fp) -> Self {
        Fp2 { c0, c1 }
    }

    /// The additive identity `0`.
    pub const ZERO: Self = Fp2 {
        c0: Fp::ZERO,
        c1: Fp::ZERO,
    };

    /// The multiplicative identity `1`.
    pub const ONE: Self = Fp2 {
        c0: Fp::ONE,
        c1: Fp::ZERO,
    };

    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    /// Double the element (`2*self`).
    pub fn double(&self) -> Self {
        Fp2 {
            c0: self.c0.double(),
            c1: self.c1.double(),
        }
    }

    /// Multiply by the non-residue `ξ = 1 + u` used to build the `Fp6`/`Fp12`
    /// tower on top of `Fp2`.
    pub fn mul_by_nonresidue(&self) -> Self {
        // (c0 + c1 u)(1 + u) = (c0 - c1) + (c0 + c1) u
        Fp2 {
            c0: &self.c0 - &self.c1,
            c1: &self.c0 + &self.c1,
        }
    }

    /// Square the element.
    pub fn square(&self) -> Self {
        // (c0 + c1 u)^2 = (c0^2 - c1^2) + 2 c0 c1 u
        //   c0' = (c0 + c1)(c0 - c1)     (= c0^2 - c1^2, since u^2 = -1)
        //   c1' = 2 c0 c1
        let c0 = &(&self.c0 + &self.c1) * &(&self.c0 - &self.c1);
        let c0c1 = &self.c0 * &self.c1;
        let c1 = &c0c1 + &c0c1;
        Fp2 { c0, c1 }
    }

    /// The multiplicative inverse.
    ///
    /// `1/(c0 + c1 u) = (c0 - c1 u) / (c0^2 + c1^2)`, where the denominator is
    /// the field norm (an [`Fp`] element). Panics if `self` is zero.
    pub fn inverse(&self) -> Self {
        let norm = &self.c0.square() + &self.c1.square();
        let inv = norm.inverse();
        Fp2 {
            c0: &self.c0 * &inv,
            c1: &(-&self.c1) * &inv,
        }
    }

    /// Sign of the element, used for point compression.
    ///
    /// Follows the usual lexicographic convention: the sign of `c0` when it is
    /// non-zero, otherwise the sign of `c1`. This flips under negation for any
    /// non-zero element, which is what point decompression relies on.
    pub fn sign(&self) -> Sign {
        if self.c0.is_zero() {
            self.c1.sign()
        } else {
            self.c0.sign()
        }
    }

    /// Raise to the power of the big-endian exponent `exp` (public exponent).
    fn power(&self, exp: &[u8]) -> Self {
        let mut a = self.clone();
        let mut q = Fp2::ONE;
        for limb in exp.iter().rev() {
            for i in 0..8 {
                if limb & (1 << i) != 0 {
                    q = &q * &a;
                }
                a = a.square();
            }
        }
        q
    }

    /// Compute a square root `x` such that `x*x == self`, if one exists.
    ///
    /// Uses the complex-method shortcut (Algorithm 9 of
    /// <https://eprint.iacr.org/2012/685.pdf>) valid because the base field has
    /// `p = 3 (mod 4)`. The returned [`CtOption`] is present only when the
    /// candidate squares back to `self`.
    pub fn sqrt(&self) -> CtOption<Self> {
        // a1 = self^((p-3)/4)
        let a1 = self.power(&P_MINUS3_DIV4_BYTES);
        // alpha = a1^2 * self = self^((p-1)/2)
        let alpha = &a1.square() * self;
        // x0 = a1 * self = self^((p+1)/4)
        let x0 = &a1 * self;

        // when alpha == -1, self is a root of unity of order p-1 and the root
        // is x0 multiplied by u (a square root of -1)
        let u = Fp2::new(Fp::zero(), Fp::one());
        let cand_neg1 = &x0 * &u;
        // otherwise the root is (1 + alpha)^((p-1)/2) * x0
        let b = (&Fp2::ONE + &alpha).power(&P_MINUS1_DIV2_BYTES);
        let cand = &b * &x0;

        let is_neg1 = alpha.ct_eq(&(-Fp2::ONE));
        let root = Fp2::ct_select(is_neg1, &cand_neg1, &cand);
        let ok = root.square().ct_eq(self);
        CtOption::from((ok, root))
    }

    /// Initialize from the `c0 || c1` big-endian bytes representation, without
    /// checking the components are canonical.
    pub const fn from_bytes_unchecked(bytes: &[u8; Self::SIZE_BYTES]) -> Self {
        let mut c0 = [0u8; Fp::SIZE_BYTES];
        let mut c1 = [0u8; Fp::SIZE_BYTES];
        let mut i = 0;
        while i < Fp::SIZE_BYTES {
            c0[i] = bytes[i];
            c1[i] = bytes[Fp::SIZE_BYTES + i];
            i += 1;
        }
        Fp2 {
            c0: Fp::from_bytes_unchecked(&c0),
            c1: Fp::from_bytes_unchecked(&c1),
        }
    }

    /// Output the `c0 || c1` big-endian bytes representation.
    pub fn to_bytes(&self) -> [u8; Self::SIZE_BYTES] {
        let mut out = [0u8; Self::SIZE_BYTES];
        out[..Fp::SIZE_BYTES].copy_from_slice(&self.c0.to_bytes());
        out[Fp::SIZE_BYTES..].copy_from_slice(&self.c1.to_bytes());
        out
    }
}

impl std::fmt::Debug for Fp2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} + {}*u", self.c0, self.c1)
    }
}

impl std::fmt::Display for Fp2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} + {}*u", self.c0, self.c1)
    }
}

impl PartialEq for Fp2 {
    fn eq(&self, other: &Self) -> bool {
        self.ct_eq(other).is_true()
    }
}
impl Eq for Fp2 {}

impl CtZero for Fp2 {
    fn ct_zero(&self) -> Choice {
        self.c0.ct_zero() & self.c1.ct_zero()
    }
    fn ct_nonzero(&self) -> Choice {
        self.c0.ct_nonzero() | self.c1.ct_nonzero()
    }
}

impl CtEqual for Fp2 {
    fn ct_eq(&self, other: &Fp2) -> Choice {
        self.c0.ct_eq(&other.c0) & self.c1.ct_eq(&other.c1)
    }
}

impl CtSelect for Fp2 {
    fn ct_select(cond: Choice, a: &Fp2, b: &Fp2) -> Fp2 {
        Fp2 {
            c0: Fp::ct_select(cond, &a.c0, &b.c0),
            c1: Fp::ct_select(cond, &a.c1, &b.c1),
        }
    }
    fn ct_assign(&mut self, cond: Choice, other: &Fp2) {
        self.c0.ct_assign(cond, &other.c0);
        self.c1.ct_assign(cond, &other.c1);
    }
}

impl From<u64> for Fp2 {
    fn from(v: u64) -> Fp2 {
        Fp2 {
            c0: Fp::from(v),
            c1: Fp::zero(),
        }
    }
}

// ****************
// Negation
// ****************

impl Neg for &Fp2 {
    type Output = Fp2;
    fn neg(self) -> Fp2 {
        Fp2 {
            c0: -&self.c0,
            c1: -&self.c1,
        }
    }
}

impl Neg for Fp2 {
    type Output = Fp2;
    fn neg(self) -> Fp2 {
        -&self
    }
}

// ****************
// Addition
// ****************

impl<'a, 'b> Add<&'b Fp2> for &'a Fp2 {
    type Output = Fp2;
    fn add(self, other: &'b Fp2) -> Fp2 {
        Fp2 {
            c0: &self.c0 + &other.c0,
            c1: &self.c1 + &other.c1,
        }
    }
}

impl<'a> Add<Fp2> for &'a Fp2 {
    type Output = Fp2;
    fn add(self, other: Fp2) -> Fp2 {
        self + &other
    }
}

impl<'b> Add<&'b Fp2> for Fp2 {
    type Output = Fp2;
    fn add(self, other: &'b Fp2) -> Fp2 {
        &self + other
    }
}

impl Add<Fp2> for Fp2 {
    type Output = Fp2;
    fn add(self, other: Fp2) -> Fp2 {
        &self + &other
    }
}

// ****************
// Subtraction
// ****************

impl<'a, 'b> Sub<&'b Fp2> for &'a Fp2 {
    type Output = Fp2;
    fn sub(self, other: &'b Fp2) -> Fp2 {
        Fp2 {
            c0: &self.c0 - &other.c0,
            c1: &self.c1 - &other.c1,
        }
    }
}

impl<'a> Sub<Fp2> for &'a Fp2 {
    type Output = Fp2;
    fn sub(self, other: Fp2) -> Fp2 {
        self - &other
    }
}

impl<'b> Sub<&'b Fp2> for Fp2 {
    type Output = Fp2;
    fn sub(self, other: &'b Fp2) -> Fp2 {
        &self - other
    }
}

impl Sub<Fp2> for Fp2 {
    type Output = Fp2;
    fn sub(self, other: Fp2) -> Fp2 {
        &self - &other
    }
}

// ****************
// Multiplication
// ****************

impl<'a, 'b> Mul<&'b Fp2> for &'a Fp2 {
    type Output = Fp2;
    fn mul(self, other: &'b Fp2) -> Fp2 {
        // (a0 + a1 u)(b0 + b1 u) = (a0 b0 - a1 b1) + (a0 b1 + a1 b0) u
        let v0 = &self.c0 * &other.c0;
        let v1 = &self.c1 * &other.c1;
        let c0 = &v0 - &v1;
        // Karatsuba for the middle term: (a0 + a1)(b0 + b1) - v0 - v1
        let c1 = &(&(&self.c0 + &self.c1) * &(&other.c0 + &other.c1)) - &(&v0 + &v1);
        Fp2 { c0, c1 }
    }
}

impl<'b> Mul<&'b Fp2> for Fp2 {
    type Output = Fp2;
    fn mul(self, other: &'b Fp2) -> Fp2 {
        &self * other
    }
}

impl<'a> Mul<Fp2> for &'a Fp2 {
    type Output = Fp2;
    fn mul(self, other: Fp2) -> Fp2 {
        self * &other
    }
}

impl Mul<Fp2> for Fp2 {
    type Output = Fp2;
    fn mul(self, other: Fp2) -> Fp2 {
        &self * &other
    }
}

impl Field for Fp2 {
    const ZERO: Fp2 = Fp2::ZERO;
    const ONE: Fp2 = Fp2::ONE;

    fn is_zero(&self) -> bool {
        self.is_zero()
    }
    fn sign(&self) -> Sign {
        self.sign()
    }
    fn double(&self) -> Fp2 {
        self.double()
    }
    fn inverse(&self) -> Fp2 {
        self.inverse()
    }
    fn square(&self) -> Fp2 {
        self.square()
    }
    fn cube(&self) -> Fp2 {
        &self.square() * self
    }
}

impl FieldSqrt for Fp2 {
    fn sqrt(&self) -> CtOption<Fp2> {
        self.sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::Fp;
    use super::Fp2;
    use crate::curve::field::Field;

    fn e(c0: u64, c1: u64) -> Fp2 {
        Fp2::new(Fp::from_u64(c0), Fp::from_u64(c1))
    }

    #[test]
    fn add_sub() {
        let a = e(3, 7);
        let b = e(10, 2);
        assert_eq!(&(&a + &b) - &b, a);
        assert_eq!(&a - &a, Fp2::ZERO);
        assert_eq!(&a + (-&a), Fp2::ZERO);
    }

    #[test]
    fn mul_matches_definition() {
        // (0 + u) * (0 + u) = u^2 = -1
        let u = e(0, 1);
        assert_eq!(&u * &u, -Fp2::ONE);
        // (a0 + a1 u)(b0 + b1 u) = (a0 b0 - a1 b1) + (a0 b1 + a1 b0) u
        let a = e(2, 3);
        let b = e(5, 7);
        // (10 - 21) + (14 + 15)u = -11 + 29u
        assert_eq!(&a * &b, Fp2::new(-Fp::from_u64(11), Fp::from_u64(29)));
    }

    #[test]
    fn square_matches_mul() {
        for (c0, c1) in [(1u64, 0), (0, 1), (2, 3), (7, 11), (123, 456)] {
            let a = e(c0, c1);
            assert_eq!(a.square(), &a * &a);
        }
    }

    #[test]
    fn inverse_roundtrip() {
        for (c0, c1) in [(1u64, 0), (0, 1), (2, 3), (7, 11), (123, 456)] {
            let a = e(c0, c1);
            assert_eq!(&a * &a.inverse(), Fp2::ONE);
        }
    }

    #[test]
    fn cube_matches() {
        let a = e(4, 9);
        assert_eq!(a.cube(), &(&a * &a) * &a);
    }

    #[test]
    fn sqrt_roundtrip() {
        // squares always have a square root, recovered up to sign
        for (c0, c1) in [(2u64, 3), (7, 11), (1, 1), (123, 456), (5, 0), (0, 8)] {
            let a = e(c0, c1);
            let sq = a.square();
            let r = sq.sqrt().into_option().expect("square has a root");
            assert_eq!(r.square(), sq, "sqrt(a^2)^2 != a^2 for ({},{})", c0, c1);
            assert!(
                r == a || r == -&a,
                "sqrt(a^2) is not +/-a for ({},{})",
                c0,
                c1
            );
        }
    }

    #[test]
    fn bytes_roundtrip() {
        let a = e(0x1234_5678, 0x9abc_def0);
        let bytes = a.to_bytes();
        assert_eq!(Fp2::from_bytes_unchecked(&bytes), a);
    }

    #[test]
    fn sign_flips_under_negation() {
        for (c0, c1) in [(3u64, 0), (0, 5), (7, 9)] {
            let a = e(c0, c1);
            assert_ne!(a.sign(), (-&a).sign());
        }
    }
}
