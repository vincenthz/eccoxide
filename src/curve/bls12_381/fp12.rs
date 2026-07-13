//! BLS12-381 degree-12 extension field `Fp12`
//!
//! `Fp12 = Fp6[w]/(w^2 - v)`, so elements are `c0 + c1*w` with `c0, c1 ∈ Fp6`
//! and `w^2 = v`. This is the codomain of the pairing (the pairing value lives
//! in the order-`r` subgroup of `Fp12*`).
//!
//! As with [`super::fp6`], only the arithmetic the pairing needs is provided
//! (add/sub/neg, mul, square, conjugate, inverse and a big-endian
//! exponentiation for the final exponentiation).

use super::fp6::Fp6;
use std::ops::{Add, Mul, Neg, Sub};

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
