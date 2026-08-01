//! BLS12-381 degree-12 extension field `Fp12`
//!
//! `Fp12 = Fp6[w]/(w^2 - v)`, so elements are `c0 + c1*w` with `c0, c1 ∈ Fp6`
//! and `w^2 = v`. This is the codomain of the pairing (the pairing value lives
//! in the order-`r` subgroup of `Fp12*`).
//!
//! As with [`super::fp6`], only the arithmetic the pairing needs is provided
//! (add/sub/neg, mul, square, conjugate, inverse and a big-endian
//! exponentiation for the final exponentiation).

use super::fp2::Fp2;
use super::fp6::Fp6;
use crate::params::bls12_381::FROBENIUS_FP12_C1_BYTES;
use std::ops::{Add, Mul, Neg, Sub};

/// `ξ^((p-1)/6)`, so that `w^p = w · FROBENIUS_C1`.
const FROBENIUS_C1: Fp2 = Fp2::from_bytes_unchecked(&FROBENIUS_FP12_C1_BYTES);

/// Element `c0 + c1*w` of `Fp12 = Fp6[w]/(w^2 - v)`.
#[derive(Clone)]
pub struct Fp12 {
    pub c0: Fp6,
    pub c1: Fp6,
}

impl Fp12 {
    pub const fn new(c0: Fp6, c1: Fp6) -> Self {
        Fp12 { c0, c1 }
    }

    pub const ZERO: Self = Fp12 {
        c0: Fp6::ZERO,
        c1: Fp6::ZERO,
    };

    pub const ONE: Self = Fp12 {
        c0: Fp6::ONE,
        c1: Fp6::ZERO,
    };

    /// Embed an `Fp6` element as `c0 + 0*w`.
    pub const fn from_fp6(c0: Fp6) -> Self {
        Fp12 { c0, c1: Fp6::ZERO }
    }

    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    /// Conjugate over `Fp6`: `c0 + c1*w -> c0 - c1*w`.
    ///
    /// This is the `p^6`-power Frobenius, and coincides with the inverse for
    /// elements of the order-`r` subgroup (they lie in the cyclotomic
    /// subgroup, where `x^{p^6} = x^{-1}`).
    pub fn conjugate(&self) -> Self {
        Fp12 {
            c0: self.c0.clone(),
            c1: -&self.c1,
        }
    }

    pub fn square(&self) -> Self {
        // (c0 + c1 w)^2 = (c0^2 + c1^2·v) + 2 c0 c1 w
        let c0c1 = &self.c0 * &self.c1;
        let c0 = &self.c0.square() + &self.c1.square().mul_by_nonresidue();
        let c1 = &c0c1 + &c0c1;
        Fp12 { c0, c1 }
    }

    /// The `p`-power Frobenius endomorphism, i.e. `self^p`.
    ///
    /// Acts coefficient-wise on the two `Fp6` components, plus the basis change
    /// `w^p = w·ξ^((p-1)/6)` coming from `w⁶ = ξ`.
    pub fn frobenius_map(&self) -> Self {
        Fp12 {
            c0: self.c0.frobenius_map(),
            c1: &self.c1.frobenius_map() * &Fp6::from_fp2(FROBENIUS_C1),
        }
    }

    /// Square an element of the cyclotomic subgroup `G_Φ6(p²)`.
    ///
    /// Granger-Scott squaring (<https://eprint.iacr.org/2009/565.pdf>): for
    /// elements satisfying the cyclotomic relation — which is the case for
    /// everything after the easy part of the final exponentiation — squaring
    /// costs 9 `Fp2` multiplications rather than the 18 of [`Self::square`].
    ///
    /// The result is **not** the square of a general `Fp12` element; only use
    /// this on the cyclotomic subgroup.
    pub fn cyclotomic_square(&self) -> Self {
        // Fp12 = Fp4 ⊕ Fp4·w ⊕ Fp4·w², where Fp4 = Fp2[s]/(s² - ξ) with s = w³.
        // The three Fp4 components are therefore the coefficient pairs
        // (c0.c0, c1.c1), (c1.c0, c0.c2) and (c0.c1, c1.c2).
        //
        // Squaring `a + b·s` in Fp4 is `(a² + ξb²) + 2ab·s`.
        fn fp4_square(a: &Fp2, b: &Fp2) -> (Fp2, Fp2) {
            let aa = a.square();
            let bb = b.square();
            let c0 = &bb.mul_by_nonresidue() + &aa;
            let c1 = &(&(a + b).square() - &aa) - &bb;
            (c0, c1)
        }

        let (t0, t1) = fp4_square(&self.c0.c0, &self.c1.c1);
        let (t2, t3) = fp4_square(&self.c1.c0, &self.c0.c2);
        let (t4, t5) = fp4_square(&self.c0.c1, &self.c1.c2);

        // each output coefficient is `3·t - 2·z` (or `3·t + 2·z`), the
        // compressed-squaring identity of the cyclotomic subgroup
        let z0 = &(&t0 - &self.c0.c0).double() + &t0;
        let z1 = &(&t1 + &self.c1.c1).double() + &t1;
        let z4 = &(&t2 - &self.c0.c1).double() + &t2;
        let z5 = &(&t3 + &self.c1.c2).double() + &t3;
        let t = t5.mul_by_nonresidue();
        let z2 = &(&t + &self.c1.c0).double() + &t;
        let z3 = &(&t4 - &self.c0.c2).double() + &t4;

        Fp12 {
            c0: Fp6::new(z0, z4, z3),
            c1: Fp6::new(z2, z1, z5),
        }
    }

    /// Multiply by the sparse element with only the `Fp2` coefficients `0`, `1`
    /// and `4` set, i.e. by `(c0 + c1*v) + (c4*v)*w`.
    ///
    /// This is the shape of the pairing line functions (see
    /// [`super::pairing`]): 13 `Fp2` multiplications instead of the 18 of a
    /// full `Fp12` multiplication.
    pub fn mul_by_014(&self, c0: &Fp2, c1: &Fp2, c4: &Fp2) -> Self {
        // Same structure as the general multiplication, with `b0 = (c0, c1, 0)`
        // and `b1 = (0, c4, 0)` handled by the sparse `Fp6` routines.
        let aa = self.c0.mul_by_01(c0, c1);
        let bb = self.c1.mul_by_1(c4);
        let o = c1 + c4;
        let c1 = &(&self.c0 + &self.c1).mul_by_01(c0, &o) - &(&aa + &bb);
        Fp12 {
            c0: &aa + &bb.mul_by_nonresidue(),
            c1,
        }
    }

    /// The multiplicative inverse (panics on zero).
    pub fn inverse(&self) -> Self {
        // (c0 + c1 w)(c0 - c1 w) = c0^2 - c1^2·v =: norm ∈ Fp6
        let norm = &self.c0.square() - &self.c1.square().mul_by_nonresidue();
        let t = norm.inverse();
        Fp12 {
            c0: &self.c0 * &t,
            c1: &(-&self.c1) * &t,
        }
    }

    /// Raise to the power of the big-endian exponent `exp` (public exponent).
    pub fn pow_bytes(&self, exp: &[u8]) -> Self {
        let mut acc = Fp12::ONE;
        let mut started = false;
        for byte in exp.iter() {
            for i in (0..8).rev() {
                if started {
                    acc = acc.square();
                }
                if (byte >> i) & 1 == 1 {
                    if started {
                        acc = &acc * self;
                    } else {
                        acc = self.clone();
                        started = true;
                    }
                }
            }
        }
        acc
    }
}

impl std::fmt::Debug for Fp12 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({:?}) + ({:?})*w", self.c0, self.c1)
    }
}

impl PartialEq for Fp12 {
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1
    }
}
impl Eq for Fp12 {}

impl Neg for &Fp12 {
    type Output = Fp12;
    fn neg(self) -> Fp12 {
        Fp12 {
            c0: -&self.c0,
            c1: -&self.c1,
        }
    }
}

impl<'a, 'b> Add<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;
    fn add(self, other: &'b Fp12) -> Fp12 {
        Fp12 {
            c0: &self.c0 + &other.c0,
            c1: &self.c1 + &other.c1,
        }
    }
}

impl<'a, 'b> Sub<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;
    fn sub(self, other: &'b Fp12) -> Fp12 {
        Fp12 {
            c0: &self.c0 - &other.c0,
            c1: &self.c1 - &other.c1,
        }
    }
}

impl<'a, 'b> Mul<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;
    fn mul(self, other: &'b Fp12) -> Fp12 {
        // (a0 + a1 w)(b0 + b1 w) = (a0 b0 + a1 b1·v) + (a0 b1 + a1 b0) w
        let aa = &self.c0 * &other.c0;
        let bb = &self.c1 * &other.c1;
        let c0 = &aa + &bb.mul_by_nonresidue();
        // (a0 + a1)(b0 + b1) - aa - bb
        let c1 = &(&(&self.c0 + &self.c1) * &(&other.c0 + &other.c1)) - &(&aa + &bb);
        Fp12 { c0, c1 }
    }
}

#[cfg(test)]
mod tests {
    use super::Fp12;
    use crate::curve::bls12_381::fp::Fp;
    use crate::curve::bls12_381::fp2::Fp2;
    use crate::curve::bls12_381::fp6::Fp6;
    use crate::params::bls12_381::P_BYTES;

    fn sample() -> Fp12 {
        let f2 = |a, b| Fp2::new(Fp::from_u64(a), Fp::from_u64(b));
        Fp12::new(
            Fp6::new(f2(1, 2), f2(3, 4), f2(5, 6)),
            Fp6::new(f2(7, 8), f2(9, 10), f2(11, 12)),
        )
    }

    #[test]
    fn mul_one() {
        let a = sample();
        assert_eq!(&a * &Fp12::ONE, a);
    }

    #[test]
    fn square_matches_mul() {
        let a = sample();
        assert_eq!(a.square(), &a * &a);
    }

    #[test]
    fn inverse_roundtrip() {
        let a = sample();
        assert_eq!(&a * &a.inverse(), Fp12::ONE);
    }

    #[test]
    fn conjugate_involution() {
        let a = sample();
        assert_eq!(a.conjugate().conjugate(), a);
    }

    #[test]
    fn frobenius_map_matches_power_p() {
        // the p-power Frobenius must agree with actual exponentiation by p
        let a = sample();
        assert_eq!(a.frobenius_map(), a.pow_bytes(&P_BYTES));
    }

    #[test]
    fn frobenius_map_six_times_is_conjugation() {
        // Frobenius^6 is the p^6-power map, which is conjugation over Fp6
        let a = sample();
        let mut f = a.clone();
        for _ in 0..6 {
            f = f.frobenius_map();
        }
        assert_eq!(f, a.conjugate());
        // and Frobenius^12 is the identity
        let mut f = a.clone();
        for _ in 0..12 {
            f = f.frobenius_map();
        }
        assert_eq!(f, a);
    }

    /// Map a sample into the cyclotomic subgroup, the way the easy part of the
    /// final exponentiation does: `a^((p^6 - 1)(p^2 + 1))`.
    fn cyclotomic_sample() -> Fp12 {
        let a = sample();
        let t = &a.conjugate() * &a.inverse();
        &t.frobenius_map().frobenius_map() * &t
    }

    #[test]
    fn cyclotomic_square_matches_square() {
        let a = cyclotomic_sample();
        // sanity: the sample really is in the cyclotomic subgroup, so its
        // conjugate is its inverse
        assert_eq!(&a * &a.conjugate(), Fp12::ONE, "sample is not cyclotomic");
        assert_eq!(a.cyclotomic_square(), a.square());
        // and it stays there under squaring
        let a2 = a.cyclotomic_square();
        assert_eq!(a2.cyclotomic_square(), a2.square());
    }

    #[test]
    fn mul_by_014_matches_mul() {
        let f2 = |a, b| Fp2::new(Fp::from_u64(a), Fp::from_u64(b));
        let a = sample();
        let (c0, c1, c4) = (f2(13, 14), f2(15, 16), f2(17, 18));
        let sparse = Fp12::new(
            Fp6::new(c0.clone(), c1.clone(), Fp2::ZERO),
            Fp6::new(Fp2::ZERO, c4.clone(), Fp2::ZERO),
        );
        assert_eq!(a.mul_by_014(&c0, &c1, &c4), &a * &sparse);
    }

    #[test]
    fn pow_bytes_matches_repeated_mul() {
        let a = sample();
        // a^13 by repeated multiplication
        let mut expected = Fp12::ONE;
        for _ in 0..13 {
            expected = &expected * &a;
        }
        assert_eq!(a.pow_bytes(&[13u8]), expected);
        // a^0 == 1, a^1 == a
        assert_eq!(a.pow_bytes(&[0u8]), Fp12::ONE);
        assert_eq!(a.pow_bytes(&[1u8]), a);
    }
}
