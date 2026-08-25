//! BLS12-381 cubic extension field `Fp6`
//!
//! `Fp6 = Fp2[v]/(v^3 - ξ)` with `ξ = 1 + u`, so elements are
//! `c0 + c1*v + c2*v^2` with `c0, c1, c2 ∈ Fp2` and `v^3 = ξ`. This is the
//! middle layer of the `Fp12` tower used by the pairing.
//!
//! Only the arithmetic needed by the pairing is provided (add/sub/neg, mul,
//! square, multiply-by-`v`, inverse); it is deliberately kept internal and does
//! not implement the [`Field`](crate::curve::field::Field) trait, since `Fp6`
//! is never used as a curve base field.

use super::fp2::Fp2;
use crate::params::bls12_381::{FROBENIUS_FP6_C1_BYTES, FROBENIUS_FP6_C2_BYTES};
use core::ops::{Add, Mul, Neg, Sub};

/// `ξ^((p-1)/3)`, so that `v^p = v · FROBENIUS_C1`.
const FROBENIUS_C1: Fp2 = Fp2::from_bytes_unchecked(&FROBENIUS_FP6_C1_BYTES);
/// `ξ^(2(p-1)/3)`, so that `(v²)^p = v² · FROBENIUS_C2`.
const FROBENIUS_C2: Fp2 = Fp2::from_bytes_unchecked(&FROBENIUS_FP6_C2_BYTES);

/// Element `c0 + c1*v + c2*v^2` of `Fp6 = Fp2[v]/(v^3 - (1 + u))`.
#[derive(Clone)]
pub struct Fp6 {
    pub c0: Fp2,
    pub c1: Fp2,
    pub c2: Fp2,
}

impl Fp6 {
    pub const fn new(c0: Fp2, c1: Fp2, c2: Fp2) -> Self {
        Fp6 { c0, c1, c2 }
    }

    pub const ZERO: Self = Fp6 {
        c0: Fp2::ZERO,
        c1: Fp2::ZERO,
        c2: Fp2::ZERO,
    };

    pub const ONE: Self = Fp6 {
        c0: Fp2::ONE,
        c1: Fp2::ZERO,
        c2: Fp2::ZERO,
    };

    /// Embed an `Fp2` element as `c0 + 0*v + 0*v^2`.
    pub const fn from_fp2(c0: Fp2) -> Self {
        Fp6 {
            c0,
            c1: Fp2::ZERO,
            c2: Fp2::ZERO,
        }
    }

    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero() && self.c2.is_zero()
    }

    /// Multiply by `v` (the non-residue of the `Fp12/Fp6` extension).
    ///
    /// `(c0 + c1 v + c2 v^2) * v = c2*v^3 + c0*v + c1*v^2 = c2*ξ + c0*v + c1*v^2`.
    pub fn mul_by_nonresidue(&self) -> Self {
        Fp6 {
            c0: self.c2.mul_by_nonresidue(),
            c1: self.c0.clone(),
            c2: self.c1.clone(),
        }
    }

    /// Square the element.
    pub fn square(&self) -> Self {
        // Chung-Hasan SQR3 (<https://eprint.iacr.org/2006/471>): three `Fp2`
        // squarings and two `Fp2` multiplications, against the six
        // multiplications a Karatsuba product with itself would cost.
        //
        // Expanding the square directly gives
        //   (a0 + a1 v + a2 v^2)^2 = (a0^2 + xi*2a1a2)
        //                          + (2a0a1 + xi*a2^2) v
        //                          + (a1^2 + 2a0a2) v^2
        // and the v^2 coefficient is recovered without its own multiplication
        // from the square of the alternating sum, since
        //   (a0 - a1 + a2)^2 = a0^2 + a1^2 + a2^2 - 2a0a1 + 2a0a2 - 2a1a2
        // so a1^2 + 2a0a2 = s2 + s1 + s3 - s0 - s4.
        let s0 = self.c0.square(); // a0^2
        let s1 = (&self.c0 * &self.c1).double(); // 2a0a1
        let s2 = (&(&self.c0 - &self.c1) + &self.c2).square(); // (a0 - a1 + a2)^2
        let s3 = (&self.c1 * &self.c2).double(); // 2a1a2
        let s4 = self.c2.square(); // a2^2

        Fp6 {
            c0: &s0 + &s3.mul_by_nonresidue(),
            c1: &s1 + &s4.mul_by_nonresidue(),
            c2: &(&(&(&s1 + &s2) + &s3) - &s0) - &s4,
        }
    }

    /// The `p`-power Frobenius endomorphism, i.e. `self^p`.
    ///
    /// Frobenius is a field automorphism fixing `Fp`, so it acts coefficient-wise
    /// on the `Fp2` components and on the basis: `v^p = v·ξ^((p-1)/3)` and
    /// `(v²)^p = v²·ξ^(2(p-1)/3)`, since `v³ = ξ`.
    pub fn frobenius_map(&self) -> Self {
        Fp6 {
            c0: self.c0.frobenius_map(),
            c1: &self.c1.frobenius_map() * &FROBENIUS_C1,
            c2: &self.c2.frobenius_map() * &FROBENIUS_C2,
        }
    }

    /// Multiply by the sparse element `c0 + c1*v` (i.e. with a zero `v²`
    /// coefficient).
    ///
    /// Five `Fp2` multiplications instead of the six of a full `Fp6`
    /// multiplication. Used by [`Fp12::mul_by_014`](super::fp12::Fp12::mul_by_014).
    pub fn mul_by_01(&self, c0: &Fp2, c1: &Fp2) -> Self {
        // Same Karatsuba as the general multiplication, with the terms that
        // involve the zero `v²` coefficient of the operand dropped:
        //   r0 = a0·c0 + ξ·(a2·c1)
        //   r1 = (a0+a1)(c0+c1) - a0·c0 - a1·c1
        //   r2 = a2·c0 + a1·c1
        let v0 = &self.c0 * c0;
        let v1 = &self.c1 * c1;

        let t = &(&(&self.c1 + &self.c2) * c1) - &v1; // a2·c1
        let r0 = &v0 + &t.mul_by_nonresidue();

        let r1 = &(&(&self.c0 + &self.c1) * &(c0 + c1)) - &(&v0 + &v1);

        let t = &(&(&self.c0 + &self.c2) * c0) - &v0; // a2·c0
        let r2 = &t + &v1;

        Fp6 {
            c0: r0,
            c1: r1,
            c2: r2,
        }
    }

    /// Multiply by the sparse element `c1*v` (only the `v` coefficient is
    /// non-zero).
    ///
    /// `(a0 + a1 v + a2 v²)·c1·v = ξ·(a2·c1) + (a0·c1)·v + (a1·c1)·v²`, i.e.
    /// three `Fp2` multiplications.
    pub fn mul_by_1(&self, c1: &Fp2) -> Self {
        Fp6 {
            c0: (&self.c2 * c1).mul_by_nonresidue(),
            c1: &self.c0 * c1,
            c2: &self.c1 * c1,
        }
    }

    /// The multiplicative inverse (panics on zero).
    pub fn inverse(&self) -> Self {
        // Standard cubic-extension inversion:
        //   t0 = c0^2 - ξ·(c1·c2)
        //   t1 = ξ·c2^2 - (c0·c1)
        //   t2 = c1^2 - (c0·c2)
        //   f  = c0·t0 + ξ·c2·t1 + ξ·c1·t2         (in Fp2)
        //   result = (t0, t1, t2) · f^{-1}
        let t0 = &self.c0.square() - &(&self.c1 * &self.c2).mul_by_nonresidue();
        let t1 = &self.c2.square().mul_by_nonresidue() - &(&self.c0 * &self.c1);
        let t2 = &self.c1.square() - &(&self.c0 * &self.c2);

        let f = &(&self.c0 * &t0) + &(&(&(&self.c2 * &t1) + &(&self.c1 * &t2))).mul_by_nonresidue();
        let f_inv = f.inverse();

        Fp6 {
            c0: &t0 * &f_inv,
            c1: &t1 * &f_inv,
            c2: &t2 * &f_inv,
        }
    }
}

impl core::fmt::Debug for Fp6 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "({}) + ({})*v + ({})*v^2", self.c0, self.c1, self.c2)
    }
}

impl PartialEq for Fp6 {
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1 && self.c2 == other.c2
    }
}
impl Eq for Fp6 {}

impl Neg for &Fp6 {
    type Output = Fp6;
    fn neg(self) -> Fp6 {
        Fp6 {
            c0: -&self.c0,
            c1: -&self.c1,
            c2: -&self.c2,
        }
    }
}

impl<'a, 'b> Add<&'b Fp6> for &'a Fp6 {
    type Output = Fp6;
    fn add(self, other: &'b Fp6) -> Fp6 {
        Fp6 {
            c0: &self.c0 + &other.c0,
            c1: &self.c1 + &other.c1,
            c2: &self.c2 + &other.c2,
        }
    }
}

impl<'a, 'b> Sub<&'b Fp6> for &'a Fp6 {
    type Output = Fp6;
    fn sub(self, other: &'b Fp6) -> Fp6 {
        Fp6 {
            c0: &self.c0 - &other.c0,
            c1: &self.c1 - &other.c1,
            c2: &self.c2 - &other.c2,
        }
    }
}

impl<'a, 'b> Mul<&'b Fp6> for &'a Fp6 {
    type Output = Fp6;
    fn mul(self, other: &'b Fp6) -> Fp6 {
        // Karatsuba over Fp2 with v^3 = ξ.
        let v0 = &self.c0 * &other.c0;
        let v1 = &self.c1 * &other.c1;
        let v2 = &self.c2 * &other.c2;

        // c0 = v0 + ξ·((c1+c2)(o1+o2) - v1 - v2)
        let t = &(&(&self.c1 + &self.c2) * &(&other.c1 + &other.c2)) - &(&v1 + &v2);
        let c0 = &v0 + &t.mul_by_nonresidue();

        // c1 = (c0+c1)(o0+o1) - v0 - v1 + ξ·v2
        let t = &(&(&self.c0 + &self.c1) * &(&other.c0 + &other.c1)) - &(&v0 + &v1);
        let c1 = &t + &v2.mul_by_nonresidue();

        // c2 = (c0+c2)(o0+o2) - v0 - v2 + v1
        let t = &(&(&self.c0 + &self.c2) * &(&other.c0 + &other.c2)) - &(&v0 + &v2);
        let c2 = &t + &v1;

        Fp6 { c0, c1, c2 }
    }
}

#[cfg(test)]
mod tests {
    use super::Fp6;
    use crate::curve::bls12_381::fp::Fp;
    use crate::curve::bls12_381::fp2::Fp2;

    fn e(a: u64, b: u64, c: u64, d: u64, f: u64, g: u64) -> Fp6 {
        Fp6::new(
            Fp2::new(Fp::from_u64(a), Fp::from_u64(b)),
            Fp2::new(Fp::from_u64(c), Fp::from_u64(d)),
            Fp2::new(Fp::from_u64(f), Fp::from_u64(g)),
        )
    }

    #[test]
    fn add_sub() {
        let a = e(1, 2, 3, 4, 5, 6);
        let b = e(7, 8, 9, 10, 11, 12);
        assert_eq!(&(&a + &b) - &b, a);
        assert_eq!(&a - &a, Fp6::ZERO);
    }

    #[test]
    fn mul_one() {
        let a = e(1, 2, 3, 4, 5, 6);
        assert_eq!(&a * &Fp6::ONE, a);
    }

    #[test]
    fn v_cubed_is_nonresidue() {
        // v * v * v == ξ (embedded in Fp6 as its constant term)
        let v = Fp6::new(Fp2::ZERO, Fp2::ONE, Fp2::ZERO);
        let xi = Fp2::ONE.mul_by_nonresidue(); // ξ = 1 + u
        assert_eq!(&(&v * &v) * &v, Fp6::from_fp2(xi));
    }

    #[test]
    fn square_matches_mul() {
        // the SQR3 formula recovers the v^2 coefficient from (a0 - a1 + a2)^2,
        // so cover the patterns where that alternating sum degenerates
        for a in [
            e(1, 2, 3, 4, 5, 6),
            Fp6::ZERO,
            Fp6::ONE,
            e(0, 0, 3, 4, 5, 6),
            e(1, 2, 0, 0, 5, 6),
            e(1, 2, 3, 4, 0, 0),
            e(1, 2, 1, 2, 1, 2), // a0 - a1 + a2 == a0
            e(1, 2, 2, 4, 1, 2), // a0 - a1 + a2 == 0
        ] {
            assert_eq!(a.square(), &a * &a, "square != mul for {:?}", a);
        }

        // and a deep chain, so any error compounds instead of cancelling
        let mut x = e(7, 11, 13, 17, 19, 23);
        let mut y = x.clone();
        for _ in 0..16 {
            x = x.square();
            y = &y * &y;
            assert_eq!(x, y);
        }
    }

    #[test]
    fn inverse_roundtrip() {
        for a in [
            e(1, 2, 3, 4, 5, 6),
            e(0, 1, 0, 0, 2, 0),
            e(9, 0, 0, 7, 0, 3),
        ] {
            assert_eq!(&a * &a.inverse(), Fp6::ONE);
        }
    }

    fn f2(a: u64, b: u64) -> Fp2 {
        Fp2::new(Fp::from_u64(a), Fp::from_u64(b))
    }

    #[test]
    fn frobenius_map_is_an_automorphism() {
        // Frobenius fixes Fp, is multiplicative, and has order 6 on Fp6
        let a = e(1, 2, 3, 4, 5, 6);
        let b = e(7, 8, 9, 10, 11, 12);
        assert_eq!(
            (&a * &b).frobenius_map(),
            &a.frobenius_map() * &b.frobenius_map()
        );
        assert_eq!(
            (&a + &b).frobenius_map(),
            &a.frobenius_map() + &b.frobenius_map()
        );
        assert_eq!(Fp6::ONE.frobenius_map(), Fp6::ONE);

        let mut f = a.clone();
        for _ in 0..6 {
            f = f.frobenius_map();
        }
        assert_eq!(f, a);
    }

    #[test]
    fn frobenius_coefficients_are_consistent() {
        // C1 = xi^((p-1)/3) and C2 = xi^(2(p-1)/3), so C2 = C1^2 and
        // C1^3 = xi^(p-1) = xi^p / xi = conj(xi) / xi
        let c1 = super::FROBENIUS_C1;
        let c2 = super::FROBENIUS_C2;
        let xi = Fp2::ONE.mul_by_nonresidue();
        assert_eq!(c1.square(), c2);
        assert_eq!(&c1.square() * &c1, &xi.frobenius_map() * &xi.inverse());
    }

    #[test]
    fn mul_by_01_matches_mul() {
        let a = e(1, 2, 3, 4, 5, 6);
        let c0 = f2(7, 8);
        let c1 = f2(9, 10);
        let sparse = Fp6::new(c0.clone(), c1.clone(), Fp2::ZERO);
        assert_eq!(a.mul_by_01(&c0, &c1), &a * &sparse);
    }

    #[test]
    fn mul_by_1_matches_mul() {
        let a = e(1, 2, 3, 4, 5, 6);
        let c1 = f2(9, 10);
        let sparse = Fp6::new(Fp2::ZERO, c1.clone(), Fp2::ZERO);
        assert_eq!(a.mul_by_1(&c1), &a * &sparse);
    }
}
