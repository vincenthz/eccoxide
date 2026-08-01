//! BLS12-381 optimal-ate pairing `e: G1 × G2 → Fp12`
//!
//! The Miller loop keeps its accumulator `T` on the sextic twist
//! `E'(Fp2): y² = x³ + 4(1 + u)` — where `G2` already lives — so all the point
//! arithmetic stays in `Fp2`. The `G1` point is instead pushed onto
//! `E'(Fp12)` by the twist map `φ(x, y) = (x·w², y·w³)`, with `w` the `Fp12`
//! tower generator (`w² = v`, `w⁶ = ξ`). Concretely:
//!
//! * `T` is held in homogeneous projective coordinates `(X : Y : Z)`, so a
//!   doubling or addition step costs a handful of `Fp2` multiplications and no
//!   inversion at all.
//! * Each step's line function, evaluated at `φ(P)`, is sparse: only three of
//!   the twelve `Fp2` coefficients are non-zero, at positions `0`, `1` and `4`.
//!   It is folded into the accumulator with
//!   [`Fp12::mul_by_014`](Fp12::mul_by_014) (13 `Fp2` muls instead of 18).
//! * The vertical-line denominators are dropped: `x_{φ(P)} - x_T = x_P·w² - x_T`
//!   lies in `Fp2[v] = Fp6`, and every element of a proper subfield of `Fp12` is
//!   killed by the final exponentiation. Working on the twist rather than
//!   untwisting `G2` likewise only changes the Miller value by a factor in
//!   `Fp4*`, which the final exponentiation kills too.
//!
//! The final exponentiation by `(p¹² − 1)/r` is split the standard way, into an
//! easy part `(p⁶ − 1)(p² + 1)` built from Frobenius and one inversion, and a
//! hard part `(p⁴ − p² + 1)/r` evaluated with an addition chain in the seed `x`
//! based on the factorization `λ₃·(p + x)·(p² + x² - 1) + 1`, where
//! `λ₃ = (x - 1)²/3`. Everything after the easy part lies in the cyclotomic
//! subgroup, so squarings use [`Fp12::cyclotomic_square`] and inversions are
//! just conjugations.

use super::fp::Fp;
use super::fp12::Fp12;
use super::fp2::Fp2;
use super::{g1, g2};
use crate::curve::{affine, projective};
use crate::params::bls12_381::HARD_EXP_LAMBDA3_BYTES;

/// BLS12-381 seed parameter `x = -0xd201000000010000`. The Miller loop runs
/// over `|x|` and the negative sign is corrected by conjugating the result.
const BLS_X: u64 = 0xd201000000010000;
const BLS_X_IS_NEGATIVE: bool = true;

/// `3b'`, three times the `b` coefficient `4(1 + u)` of the G2 twist.
const B3_TWIST: Fp2 = Fp2::from_bytes_unchecked(&crate::params::bls12_381::g2::B3_BYTES);

/// The three non-zero `Fp2` coefficients `(c0, c1, c4)` of a line function,
/// before the `G1` coordinates are folded in (see [`ell`]).
type LineCoeffs = (Fp2, Fp2, Fp2);

/// The Miller loop accumulator: a point of the twist `E'(Fp2)` in homogeneous
/// projective coordinates `(X : Y : Z)`, i.e. the affine point `(X/Z, Y/Z)`.
///
/// This is the same representation the `G2` group uses, so the shared
/// [`projective::Point`] serves as the carrier; only the doubling and addition
/// *formulas* are specific to the pairing, since each one has to hand back the
/// line coefficients that fall out of its intermediate values.
type Point = projective::Point<Fp2>;

/// Replace `t` by `2t` and return the coefficients of the tangent line at `t`.
///
/// Homogeneous projective doubling with line evaluation (algorithm 26 of
/// <https://eprint.iacr.org/2010/354.pdf>), scaled by 4 so that no division by
/// two is needed. Both the point and the line are only defined up to a non-zero
/// `Fp2` factor — the point because it is projective, the line because such a
/// factor does not survive the final exponentiation.
fn doubling_step(t: &mut Point) -> LineCoeffs {
    let xy = &t.x * &t.y;
    let xx = t.x.square();
    let yy = t.y.square();
    let zz = t.z.square();

    let e = &B3_TWIST * &zz; // 3b'·Z²
    let f = &e.double() + &e; // 9b'·Z²
    let two_yz = &(&t.y + &t.z).square() - &(&yy + &zz); // 2YZ
    let e4 = e.square().double().double(); // 4·(3b'Z²)²

    // 2T, uniformly scaled by 4 so the halving of the usual formulas drops out
    t.x = &xy.double() * &(&yy - &f);
    t.y = &(&yy + &f).square() - &(&e4.double() + &e4);
    t.z = &yy.double().double() * &two_yz;

    // ℓ = (3b'Z² - Y²) + (3X²·x_P)·v - (2YZ·y_P)·v·w
    (&e - &yy, &xx.double() + &xx, -&two_yz)
}

/// Replace `t` by `t + q` and return the coefficients of the line through `t`
/// and `q`, with `q = (qx, qy)` affine.
fn addition_step(t: &mut Point, qx: &Fp2, qy: &Fp2) -> LineCoeffs {
    // θ and λ are the numerator and denominator of the chord's slope, scaled by
    // Z: θ/λ = (y_T - qy) / (x_T - qx).
    let theta = &t.y - &(qy * &t.z);
    let lambda = &t.x - &(qx * &t.z);
    let ll = lambda.square();
    let lll = &lambda * &ll;
    let x_ll = &t.x * &ll;
    let h = &(&lll + &(&t.z * &theta.square())) - &x_ll.double();

    // T + Q
    t.x = &lambda * &h;
    t.y = &(&theta * &(&x_ll - &h)) - &(&lll * &t.y);
    t.z = &t.z * &lll;

    // ℓ = (θ·qx - λ·qy) - (θ·x_P)·v + (λ·y_P)·v·w
    (&(&theta * qx) - &(&lambda * qy), -&theta, lambda)
}

/// Fold a line value into the Miller accumulator: `f * ℓ(φ(P))`.
///
/// The line is `c0 + (c1·x_P)·v + (c4·y_P)·v·w` once the `G1` coordinates of
/// `P` are multiplied in, which is exactly the sparse shape
/// [`Fp12::mul_by_014`] expects.
fn ell(f: &Fp12, coeffs: &LineCoeffs, px: &Fp, py: &Fp) -> Fp12 {
    f.mul_by_014(&coeffs.0, &coeffs.1.mul_by_fp(px), &coeffs.2.mul_by_fp(py))
}

/// The Miller loop `f_{|x|, Q}(φ(P))`, for `P = (px, py)` on `G1` and
/// `Q = (qx, qy)` on `G2`.
fn miller(px: &Fp, py: &Fp, qx: &Fp2, qy: &Fp2) -> Fp12 {
    let mut f = Fp12::ONE;
    let mut t = Point::from_affine(&affine::Point {
        x: qx.clone(),
        y: qy.clone(),
    });

    // |x| has bit length 64 (bit 63 set): the top bit is consumed by the
    // initial T = Q, then process bits 62..=0.
    for i in (0..63).rev() {
        let coeffs = doubling_step(&mut t);
        f = ell(&f.square(), &coeffs, px, py);

        if (BLS_X >> i) & 1 == 1 {
            let coeffs = addition_step(&mut t, qx, qy);
            f = ell(&f, &coeffs, px, py);
        }
    }

    f
}

/// Raise to the power of the big-endian public exponent `exp`, squaring with
/// [`Fp12::cyclotomic_square`].
///
/// `self` must lie in the cyclotomic subgroup, which every input here does: it
/// is closed under multiplication and everything is downstream of [`easy_part`].
fn cyclotomic_pow(f: &Fp12, exp: &[u8]) -> Fp12 {
    let mut acc = Fp12::ONE;
    let mut started = false;
    for byte in exp.iter() {
        for i in (0..8).rev() {
            if started {
                acc = acc.cyclotomic_square();
            }
            if (byte >> i) & 1 == 1 {
                if started {
                    acc = &acc * f;
                } else {
                    acc = f.clone();
                    started = true;
                }
            }
        }
    }
    acc
}

/// `f^x` for the seed `x`, which has only 6 set bits.
///
/// `x` is negative, and inversion in the cyclotomic subgroup is conjugation, so
/// the sign costs nothing.
fn exp_by_x(f: &Fp12) -> Fp12 {
    let t = cyclotomic_pow(f, &BLS_X.to_be_bytes());
    if BLS_X_IS_NEGATIVE {
        t.conjugate()
    } else {
        t
    }
}

/// The easy part of the final exponentiation: `f^((p⁶ - 1)(p² + 1))`.
///
/// `f^(p⁶ - 1)` is `conjugate(f) · f⁻¹` — the one and only inversion in the
/// whole pairing — and `f^(p² + 1)` is a double Frobenius and a multiplication.
/// The result satisfies `t^(p⁶ + 1) = 1`, i.e. it lies in the cyclotomic
/// subgroup, where `conjugate` *is* the inverse.
fn easy_part(f: &Fp12) -> Fp12 {
    let t = &f.conjugate() * &f.inverse();
    &t.frobenius_map().frobenius_map() * &t
}

/// The hard part of the final exponentiation: `m^((p⁴ - p² + 1)/r)`.
///
/// Writing the exponent in base `p` and recognizing the digits as multiples of
/// `λ₃ = (x - 1)²/3` gives the factorization
///
/// ```text
/// (p⁴ - p² + 1)/r = λ₃·(p + x)·(p² + x² - 1) + 1
/// ```
///
/// which needs one exponentiation by `λ₃` (126 bits) and three by the seed `x`
/// (64 bits, 6 set bits each), plus four Frobenius applications. Note this
/// computes the exponent *exactly*; the cheaper `(x-1)²`-based chains found in
/// several libraries evaluate three times the hard part, which is still a valid
/// pairing but cubes the value.
fn hard_part(m: &Fp12) -> Fp12 {
    // y = m^λ₃
    let y = cyclotomic_pow(m, &HARD_EXP_LAMBDA3_BYTES);
    // y1 = y^(p + x)
    let y1 = &y.frobenius_map() * &exp_by_x(&y);
    // y2 = y1^(p² + x² - 1)
    let y2 = &(&y1.frobenius_map().frobenius_map() * &exp_by_x(&exp_by_x(&y1))) * &y1.conjugate();

    m * &y2
}

/// Compute the optimal-ate pairing `e(p, q)` of a `G1` point and a `G2` point.
///
/// The inputs are affine points (i.e. not the point at infinity); the result
/// lives in the order-`r` subgroup of `Fp12*`.
pub fn pairing(p: &g1::PointAffine, q: &g2::PointAffine) -> Fp12 {
    let (px, py) = p.to_coordinate();
    let (qx, qy) = q.to_coordinate();

    let f = hard_part(&easy_part(&miller(px, py, qx, qy)));

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
    use super::{addition_step, doubling_step, pairing, Point, BLS_X};
    use crate::curve::affine;
    use crate::curve::bls12_381::fp12::Fp12;
    use crate::curve::bls12_381::fp2::Fp2;
    use crate::curve::bls12_381::scalar::Scalar;
    use crate::curve::bls12_381::{g1, g2};
    use crate::params::bls12_381::{FINAL_EXP_BYTES, ORDER_BYTES};

    fn g1_mul(a: &Scalar) -> g1::PointAffine {
        g1::Point::mul_base(a).to_affine().unwrap()
    }
    fn g2_mul(a: &Scalar) -> g2::PointAffine {
        g2::Point::mul_base(a).to_affine().unwrap()
    }

    /// The straightforward, correctness-first pairing the optimized one
    /// replaced: `G2` is untwisted onto `E(Fp12)` with `ψ(x, y) = (x/w², y/w³)`,
    /// `G1` is embedded directly, and the Miller loop is a textbook affine
    /// double-and-add over `E(Fp12)` keeping the vertical-line denominators.
    /// Kept here as an independent reference for the fast implementation.
    mod reference {
        use crate::curve::bls12_381::fp::Fp;
        use crate::curve::bls12_381::fp12::Fp12;
        use crate::curve::bls12_381::fp2::Fp2;
        use crate::curve::bls12_381::fp6::Fp6;
        use crate::curve::bls12_381::{g1, g2};
        use crate::params::bls12_381::FINAL_EXP_BYTES;

        fn fp_to_fp12(a: &Fp) -> Fp12 {
            Fp12::from_fp6(Fp6::from_fp2(Fp2::new(a.clone(), Fp::zero())))
        }

        fn fp2_to_fp12(a: &Fp2) -> Fp12 {
            Fp12::from_fp6(Fp6::from_fp2(a.clone()))
        }

        fn miller(px: &Fp12, py: &Fp12, qx: &Fp12, qy: &Fp12) -> Fp12 {
            let mut f_num = Fp12::ONE;
            let mut f_den = Fp12::ONE;
            let mut tx = qx.clone();
            let mut ty = qy.clone();

            for i in (0..63).rev() {
                // doubling: tangent at T (λ = 3·x_T² / 2·y_T, the curve has a = 0)
                let x2 = tx.square();
                let num = &(&x2 + &x2) + &x2;
                let den = &ty + &ty;
                let lambda = &num * &den.inverse();
                let x3 = &lambda.square() - &(&tx + &tx);
                let y3 = &(&lambda * &(&tx - &x3)) - &ty;
                let line = &(py - &ty) - &(&lambda * &(px - &tx));
                let vert = px - &x3;
                f_num = &f_num.square() * &line;
                f_den = &f_den.square() * &vert;
                tx = x3;
                ty = y3;

                if (super::BLS_X >> i) & 1 == 1 {
                    // addition: line through T and Q
                    let num = &ty - qy;
                    let den = &tx - qx;
                    let lambda = &num * &den.inverse();
                    let x3 = &(&lambda.square() - &tx) - qx;
                    let y3 = &(&lambda * &(&tx - &x3)) - &ty;
                    let line = &(py - &ty) - &(&lambda * &(px - &tx));
                    let vert = px - &x3;
                    f_num = &f_num * &line;
                    f_den = &f_den * &vert;
                    tx = x3;
                    ty = y3;
                }
            }

            &f_num * &f_den.inverse()
        }

        pub fn pairing(p: &g1::PointAffine, q: &g2::PointAffine) -> Fp12 {
            let (px, py) = p.to_coordinate();
            let (qx, qy) = q.to_coordinate();

            let w = Fp12::new(Fp6::ZERO, Fp6::ONE);
            let w2 = w.square();
            let w3 = &w2 * &w;

            let m = miller(
                &fp_to_fp12(px),
                &fp_to_fp12(py),
                &(&fp2_to_fp12(qx) * &w2.inverse()),
                &(&fp2_to_fp12(qy) * &w3.inverse()),
            );
            m.pow_bytes(&FINAL_EXP_BYTES).conjugate()
        }
    }

    #[test]
    fn matches_reference_implementation() {
        // the projective / sparse-line Miller loop must agree with the naive
        // affine one over E(Fp12), not merely be *some* valid pairing: the
        // subfield factors it drops have to be exactly that.
        let cases = [
            (1u64, 1u64),
            (2, 3),
            (7, 1),
            (0x9e3779b97f4a7c15, 0xdeadbeef),
        ];
        for (a, b) in cases {
            let p = g1_mul(&Scalar::from(a));
            let q = g2_mul(&Scalar::from(b));
            assert_eq!(
                pairing(&p, &q),
                reference::pairing(&p, &q),
                "e([{}]G1, [{}]G2) differs from the reference",
                a,
                b
            );
        }
    }

    #[test]
    fn final_exponentiation_matches_full_exponent() {
        // the easy+hard addition chain must compute exactly (p^12 - 1)/r, not
        // merely some multiple of it that would still be a valid pairing
        let q = g2::PointAffine::GENERATOR;
        for a in [1u64, 2, 3, 0x9e3779b97f4a7c15] {
            let p = g1_mul(&Scalar::from(a));
            let (px, py) = p.to_coordinate();
            let (qx, qy) = q.to_coordinate();
            let m = super::miller(px, py, qx, qy);
            assert_eq!(
                super::hard_part(&super::easy_part(&m)),
                m.pow_bytes(&FINAL_EXP_BYTES),
                "final exponentiation differs from f^((p^12-1)/r) for a = {}",
                a
            );
        }
    }

    #[test]
    fn miller_accumulator_tracks_scalar_mul() {
        // the projective doubling / addition steps must leave T = [|x|]Q
        let q = g2::PointAffine::GENERATOR;
        let (qx, qy) = q.to_coordinate();
        let mut t = Point {
            x: qx.clone(),
            y: qy.clone(),
            z: Fp2::ONE,
        };
        for i in (0..63).rev() {
            doubling_step(&mut t);
            if (BLS_X >> i) & 1 == 1 {
                addition_step(&mut t, qx, qy);
            }
        }

        let expected = super::Point::from_affine(&affine::Point {
            x: qx.clone(),
            y: qy.clone(),
        })
        .scale_a0::<g2::Curve>(&BLS_X.to_be_bytes());
        assert_eq!(t, expected);
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
