//! BLS12-381 optimal-ate pairing `e: G1 × G2 → Fp12`
//!
//! This is a deliberately straightforward, correctness-first implementation:
//!
//! * `G2` points are moved onto `E(Fp12)` with the untwist map
//!   `ψ(x, y) = (x/w², y/w³)` (where `w` is the `Fp12` tower generator, so
//!   `w⁶ = ξ`), and `G1` points are embedded directly. The whole Miller loop
//!   then runs as a textbook double-and-add over `E(Fp12)` with the generic
//!   Weierstrass line functions, so no specialized twist / sparse-multiplication
//!   formulas are needed.
//! * The final exponentiation is a single square-and-multiply by the full
//!   exponent `(p¹² − 1)/r`, so no Frobenius-coefficient constants are needed.
//!
//! Both choices trade speed for simplicity and auditability; the optimized
//! projective line functions, sparse `Fp12` multiplication and the
//! Frobenius-based final exponentiation are natural future optimizations.

use super::fp::Fp;
use super::fp12::Fp12;
use super::fp2::Fp2;
use super::fp6::Fp6;
use super::{g1, g2};
use crate::params::bls12_381::FINAL_EXP_BYTES;

/// BLS12-381 seed parameter `x = -0xd201000000010000`. The Miller loop runs
/// over `|x|` and the negative sign is corrected by conjugating the result.
const BLS_X: u64 = 0xd201000000010000;
const BLS_X_IS_NEGATIVE: bool = true;

/// An affine point on `E(Fp12): y² = x³ + 4`, used inside the Miller loop.
struct Fp12Point {
    x: Fp12,
    y: Fp12,
}

/// Embed an `Fp` element into `Fp12`.
fn fp_to_fp12(a: &Fp) -> Fp12 {
    Fp12::from_fp6(Fp6::from_fp2(Fp2::new(a.clone(), Fp::zero())))
}

/// Embed an `Fp2` element into `Fp12`.
fn fp2_to_fp12(a: &Fp2) -> Fp12 {
    Fp12::from_fp6(Fp6::from_fp2(a.clone()))
}

/// The Miller loop `f_{|x|, Q}(P)`, with `P` and `Q` already mapped onto
/// `E(Fp12)`. Numerator (line) and denominator (vertical) factors are
/// accumulated separately and combined with a single inversion at the end.
fn miller(p: &Fp12Point, q: &Fp12Point) -> Fp12 {
    let mut f_num = Fp12::ONE;
    let mut f_den = Fp12::ONE;
    let mut tx = q.x.clone();
    let mut ty = q.y.clone();

    // |x| has bit length 64 (bit 63 set): the top bit is consumed by the
    // initial T = Q, then process bits 62..=0.
    for i in (0..63).rev() {
        // --- doubling: tangent at T, then T = 2T ---
        // λ = 3·x_T² / (2·y_T)   (curve has a = 0)
        let x2 = tx.square();
        let num = &(&x2 + &x2) + &x2;
        let den = &ty + &ty;
        let lambda = &num * &den.inverse();
        let x3 = &lambda.square() - &(&tx + &tx);
        let y3 = &(&lambda * &(&tx - &x3)) - &ty;
        // line ℓ(P) through T (slope λ), and vertical v(P) at 2T
        let line = &(&p.y - &ty) - &(&lambda * &(&p.x - &tx));
        let vert = &p.x - &x3;
        let fn2 = f_num.square();
        let fd2 = f_den.square();
        f_num = &fn2 * &line;
        f_den = &fd2 * &vert;
        tx = x3;
        ty = y3;

        if (BLS_X >> i) & 1 == 1 {
            // --- addition: line through T and Q, then T = T + Q ---
            let num = &ty - &q.y;
            let den = &tx - &q.x;
            let lambda = &num * &den.inverse();
            let x3 = &(&lambda.square() - &tx) - &q.x;
            let y3 = &(&lambda * &(&tx - &x3)) - &ty;
            let line = &(&p.y - &ty) - &(&lambda * &(&p.x - &tx));
            let vert = &p.x - &x3;
            f_num = &f_num * &line;
            f_den = &f_den * &vert;
            tx = x3;
            ty = y3;
        }
    }

    &f_num * &f_den.inverse()
}

/// Compute the optimal-ate pairing `e(p, q)` of a `G1` point and a `G2` point.
///
/// The inputs are affine points (i.e. not the point at infinity); the result
/// lives in the order-`r` subgroup of `Fp12*`.
pub fn pairing(p: &g1::PointAffine, q: &g2::PointAffine) -> Fp12 {
    let (px, py) = p.to_coordinate();
    let (qx, qy) = q.to_coordinate();

    // embed P directly into E(Fp12)
    let p12 = Fp12Point {
        x: fp_to_fp12(px),
        y: fp_to_fp12(py),
    };

    // untwist Q onto E(Fp12): ψ(x, y) = (x/w², y/w³), where w is the tower
    // generator (w² = v, w⁶ = ξ), so this lands on y² = x³ + 4.
    let w = Fp12::new(Fp6::ZERO, Fp6::ONE);
    let w2 = w.square();
    let w3 = &w2 * &w;
    let q12 = Fp12Point {
        x: &fp2_to_fp12(qx) * &w2.inverse(),
        y: &fp2_to_fp12(qy) * &w3.inverse(),
    };

    let m = miller(&p12, &q12);
    let f = m.pow_bytes(&FINAL_EXP_BYTES);

    // x is negative: the pairing for -|x| is the inverse, which in the
    // order-r subgroup is the conjugation over Fp6.
    if BLS_X_IS_NEGATIVE {
        f.conjugate()
    } else {
        f
    }
}

#[cfg(test)]
mod tests {
    use super::pairing;
    use crate::curve::bls12_381::fp12::Fp12;
    use crate::curve::bls12_381::scalar::Scalar;
    use crate::curve::bls12_381::{g1, g2};
    use crate::params::bls12_381::ORDER_BYTES;

    fn g1_mul(a: &Scalar) -> g1::PointAffine {
        g1::Point::mul_base(a).to_affine().unwrap()
    }
    fn g2_mul(a: &Scalar) -> g2::PointAffine {
        g2::Point::mul_base(a).to_affine().unwrap()
    }

    #[test]
    fn bilinearity_and_non_degeneracy() {
        let g1 = g1::PointAffine::GENERATOR;
        let g2 = g2::PointAffine::GENERATOR;

        let base = pairing(&g1, &g2);

        // non-degenerate: e(G1, G2) != 1
        assert_ne!(base, Fp12::ONE, "pairing is degenerate");

        // lands in the order-r subgroup: e(G1, G2)^r == 1
        assert_eq!(base.pow_bytes(&ORDER_BYTES), Fp12::ONE, "e(G1,G2)^r != 1");

        let a = Scalar::from(0x9e3779b97f4a7c15u64);
        let b = Scalar::from(0x1234_5678_9abc_def0u64);

        let pa = g1_mul(&a);
        let qb = g2_mul(&b);
        let qa = g2_mul(&a);

        // linearity in each argument: e([a]G1, G2) == e(G1, [a]G2) == e(G1,G2)^a
        let base_a = base.pow_bytes(&a.to_bytes());
        assert_eq!(pairing(&pa, &g2), base_a, "e([a]G1,G2) != e(G1,G2)^a");
        assert_eq!(pairing(&g1, &qa), base_a, "e(G1,[a]G2) != e(G1,G2)^a");

        // full bilinearity: e([a]G1, [b]G2) == e(G1,G2)^(a*b)
        let ab = &a * &b;
        assert_eq!(
            pairing(&pa, &qb),
            base.pow_bytes(&ab.to_bytes()),
            "e([a]G1,[b]G2) != e(G1,G2)^(ab)"
        );
    }

    #[test]
    fn additive_in_first_argument() {
        // e(P1 + P2, Q) == e(P1, Q) * e(P2, Q)
        let g2 = g2::PointAffine::GENERATOR;
        let p1 = g1_mul(&Scalar::from(3));
        let p2 = g1_mul(&Scalar::from(5));
        let p12 = g1_mul(&Scalar::from(8));

        let lhs = pairing(&p12, &g2);
        let rhs = &pairing(&p1, &g2) * &pairing(&p2, &g2);
        assert_eq!(lhs, rhs);
    }
}
