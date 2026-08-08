//! Divan micro-benchmarks for BLS12-381.
//!
//! Covers the base field `Fp`, the `Fp2`/`Fp6`/`Fp12` extension tower — including
//! the sparse multiplications, the Frobenius maps and the cyclotomic squaring the
//! pairing is built on — the scalar field, the G1 and G2 groups, and the
//! optimal-ate pairing itself.
//!
//! Run with:
//!
//! ```text
//! cargo bench --bench bls12_381
//! cargo bench --bench bls12_381 -- pairing        # only the pairing
//! cargo bench --bench bls12_381 -- fp12           # only the Fp12 arithmetic
//! cargo bench --bench bls12_381 -- g1::scalar_mul
//! ```

fn main() {
    divan::main();
}

#[cfg(feature = "bls12-381")]
mod bls12_381 {
    // NOTE: `g1` and `g2` are deliberately not imported here — the group
    // benchmarks below generate modules of those names, which divan uses to
    // label them.
    use eccoxide::curve::bls12_381::{Fp, Fp12, Fp2, Fp6, Scalar};

    // --- sample values ----------------------------------------------------

    /// A non-trivial, full-width `Fp` element (the G1 generator's x).
    fn fp_a() -> Fp {
        eccoxide::curve::bls12_381::g1::Curve::generator().0.clone()
    }

    /// A second, distinct full-width `Fp` element (the G1 generator's y).
    fn fp_b() -> Fp {
        eccoxide::curve::bls12_381::g1::Curve::generator().1.clone()
    }

    /// A non-trivial, full-width `Fp2` element (the G2 generator's x).
    fn fp2_a() -> Fp2 {
        eccoxide::curve::bls12_381::g2::Curve::generator().0.clone()
    }

    /// A second, distinct full-width `Fp2` element (the G2 generator's y).
    fn fp2_b() -> Fp2 {
        eccoxide::curve::bls12_381::g2::Curve::generator().1.clone()
    }

    fn fp6_a() -> Fp6 {
        Fp6::new(fp2_a(), fp2_b(), &fp2_a() * &fp2_b())
    }

    fn fp6_b() -> Fp6 {
        Fp6::new(&fp2_a() + &fp2_b(), fp2_a().square(), fp2_b().square())
    }

    fn fp12_a() -> Fp12 {
        Fp12::new(fp6_a(), fp6_b())
    }

    fn fp12_b() -> Fp12 {
        Fp12::new(fp6_b(), fp6_a())
    }

    /// A non-trivial scalar (a reduced product, so always `< order`).
    fn sc_a() -> Scalar {
        let mut s = Scalar::from_u64(0x0123_4567_89ab_cdef);
        for _ in 0..5 {
            s = s.square();
        }
        s
    }

    /// A second, distinct non-trivial scalar.
    fn sc_b() -> Scalar {
        let mut s = Scalar::from_u64(0xfedc_ba98_7654_3210);
        for _ in 0..5 {
            s = s.square();
        }
        s
    }

    // --- base prime field Fp ----------------------------------------------
    mod fp {
        use super::*;
        use divan::{black_box, Bencher};

        #[divan::bench]
        fn add(bencher: Bencher) {
            let (a, b) = (fp_a(), fp_b());
            bencher.bench(|| black_box(&a) + black_box(&b));
        }

        #[divan::bench]
        fn sub(bencher: Bencher) {
            let (a, b) = (fp_a(), fp_b());
            bencher.bench(|| black_box(&a) - black_box(&b));
        }

        #[divan::bench]
        fn mul(bencher: Bencher) {
            let (a, b) = (fp_a(), fp_b());
            bencher.bench(|| black_box(&a) * black_box(&b));
        }

        #[divan::bench]
        fn square(bencher: Bencher) {
            let a = fp_a();
            bencher.bench(|| black_box(&a).square());
        }

        #[divan::bench]
        fn double(bencher: Bencher) {
            let a = fp_a();
            bencher.bench(|| black_box(&a).double());
        }

        #[divan::bench]
        fn neg(bencher: Bencher) {
            let a = fp_a();
            bencher.bench(|| -black_box(&a));
        }

        #[divan::bench]
        fn inverse(bencher: Bencher) {
            let a = fp_a();
            bencher.bench(|| black_box(&a).inverse());
        }

        #[divan::bench]
        fn sqrt(bencher: Bencher) {
            // square first so the input is guaranteed to be a residue
            let a = fp_a().square();
            bencher.bench(|| black_box(&a).sqrt());
        }

        #[divan::bench]
        fn to_bytes(bencher: Bencher) {
            let a = fp_a();
            bencher.bench(|| black_box(&a).to_bytes());
        }

        #[divan::bench]
        fn from_bytes(bencher: Bencher) {
            let bytes = fp_a().to_bytes();
            bencher.bench(|| Fp::from_bytes(black_box(&bytes)));
        }
    }

    // --- quadratic extension Fp2 (the G2 base field) ----------------------
    mod fp2 {
        use super::*;
        use divan::{black_box, Bencher};

        #[divan::bench]
        fn add(bencher: Bencher) {
            let (a, b) = (fp2_a(), fp2_b());
            bencher.bench(|| black_box(&a) + black_box(&b));
        }

        #[divan::bench]
        fn sub(bencher: Bencher) {
            let (a, b) = (fp2_a(), fp2_b());
            bencher.bench(|| black_box(&a) - black_box(&b));
        }

        #[divan::bench]
        fn mul(bencher: Bencher) {
            let (a, b) = (fp2_a(), fp2_b());
            bencher.bench(|| black_box(&a) * black_box(&b));
        }

        #[divan::bench]
        fn square(bencher: Bencher) {
            let a = fp2_a();
            bencher.bench(|| black_box(&a).square());
        }

        #[divan::bench]
        fn double(bencher: Bencher) {
            let a = fp2_a();
            bencher.bench(|| black_box(&a).double());
        }

        /// Scaling by a base field element, as the pairing line functions do.
        #[divan::bench]
        fn mul_by_fp(bencher: Bencher) {
            let (a, b) = (fp2_a(), fp_a());
            bencher.bench(|| black_box(&a).mul_by_fp(black_box(&b)));
        }

        /// Multiplication by the tower non-residue `ξ = 1 + u`.
        #[divan::bench]
        fn mul_by_nonresidue(bencher: Bencher) {
            let a = fp2_a();
            bencher.bench(|| black_box(&a).mul_by_nonresidue());
        }

        #[divan::bench]
        fn frobenius_map(bencher: Bencher) {
            let a = fp2_a();
            bencher.bench(|| black_box(&a).frobenius_map());
        }

        #[divan::bench]
        fn inverse(bencher: Bencher) {
            let a = fp2_a();
            bencher.bench(|| black_box(&a).inverse());
        }

        #[divan::bench]
        fn sqrt(bencher: Bencher) {
            let a = fp2_a().square();
            bencher.bench(|| black_box(&a).sqrt());
        }
    }

    // --- cubic extension Fp6 ----------------------------------------------
    mod fp6 {
        use super::*;
        use divan::{black_box, Bencher};

        #[divan::bench]
        fn add(bencher: Bencher) {
            let (a, b) = (fp6_a(), fp6_b());
            bencher.bench(|| black_box(&a) + black_box(&b));
        }

        #[divan::bench]
        fn mul(bencher: Bencher) {
            let (a, b) = (fp6_a(), fp6_b());
            bencher.bench(|| black_box(&a) * black_box(&b));
        }

        #[divan::bench]
        fn square(bencher: Bencher) {
            let a = fp6_a();
            bencher.bench(|| black_box(&a).square());
        }

        /// Sparse multiplication by `c0 + c1*v`, used by `Fp12::mul_by_014`.
        #[divan::bench]
        fn mul_by_01(bencher: Bencher) {
            let (a, c0, c1) = (fp6_a(), fp2_a(), fp2_b());
            bencher.bench(|| black_box(&a).mul_by_01(black_box(&c0), black_box(&c1)));
        }

        /// Sparse multiplication by `c1*v`, used by `Fp12::mul_by_014`.
        #[divan::bench]
        fn mul_by_1(bencher: Bencher) {
            let (a, c1) = (fp6_a(), fp2_b());
            bencher.bench(|| black_box(&a).mul_by_1(black_box(&c1)));
        }

        #[divan::bench]
        fn mul_by_nonresidue(bencher: Bencher) {
            let a = fp6_a();
            bencher.bench(|| black_box(&a).mul_by_nonresidue());
        }

        #[divan::bench]
        fn frobenius_map(bencher: Bencher) {
            let a = fp6_a();
            bencher.bench(|| black_box(&a).frobenius_map());
        }

        #[divan::bench]
        fn inverse(bencher: Bencher) {
            let a = fp6_a();
            bencher.bench(|| black_box(&a).inverse());
        }
    }

    // --- degree-12 extension Fp12 (the pairing codomain) ------------------
    mod fp12 {
        use super::*;
        use divan::{black_box, Bencher};

        /// An element of the cyclotomic subgroup, as the easy part of the final
        /// exponentiation produces; `cyclotomic_square` is only valid there.
        fn cyclotomic() -> Fp12 {
            let a = fp12_a();
            let t = &a.conjugate() * &a.inverse();
            &t.frobenius_map().frobenius_map() * &t
        }

        #[divan::bench]
        fn add(bencher: Bencher) {
            let (a, b) = (fp12_a(), fp12_b());
            bencher.bench(|| black_box(&a) + black_box(&b));
        }

        #[divan::bench]
        fn mul(bencher: Bencher) {
            let (a, b) = (fp12_a(), fp12_b());
            bencher.bench(|| black_box(&a) * black_box(&b));
        }

        #[divan::bench]
        fn square(bencher: Bencher) {
            let a = fp12_a();
            bencher.bench(|| black_box(&a).square());
        }

        /// Granger-Scott squaring, the workhorse of the final exponentiation.
        #[divan::bench]
        fn cyclotomic_square(bencher: Bencher) {
            let a = cyclotomic();
            bencher.bench(|| black_box(&a).cyclotomic_square());
        }

        /// The sparse multiplication the Miller loop folds each line into.
        #[divan::bench]
        fn mul_by_014(bencher: Bencher) {
            let a = fp12_a();
            let (c0, c1, c4) = (fp2_a(), fp2_b(), &fp2_a() + &fp2_b());
            bencher
                .bench(|| black_box(&a).mul_by_014(black_box(&c0), black_box(&c1), black_box(&c4)));
        }

        #[divan::bench]
        fn conjugate(bencher: Bencher) {
            let a = fp12_a();
            bencher.bench(|| black_box(&a).conjugate());
        }

        #[divan::bench]
        fn frobenius_map(bencher: Bencher) {
            let a = fp12_a();
            bencher.bench(|| black_box(&a).frobenius_map());
        }

        #[divan::bench]
        fn inverse(bencher: Bencher) {
            let a = fp12_a();
            bencher.bench(|| black_box(&a).inverse());
        }
    }

    // --- scalar field F(r) ------------------------------------------------
    mod scalar {
        use super::*;
        use divan::{black_box, Bencher};

        #[divan::bench]
        fn add(bencher: Bencher) {
            let (a, b) = (sc_a(), sc_b());
            bencher.bench(|| black_box(&a) + black_box(&b));
        }

        #[divan::bench]
        fn sub(bencher: Bencher) {
            let (a, b) = (sc_a(), sc_b());
            bencher.bench(|| black_box(&a) - black_box(&b));
        }

        #[divan::bench]
        fn mul(bencher: Bencher) {
            let (a, b) = (sc_a(), sc_b());
            bencher.bench(|| black_box(&a) * black_box(&b));
        }

        #[divan::bench]
        fn square(bencher: Bencher) {
            let a = sc_a();
            bencher.bench(|| black_box(&a).square());
        }

        #[divan::bench]
        fn inverse(bencher: Bencher) {
            let a = sc_a();
            bencher.bench(|| black_box(&a).inverse());
        }

        #[divan::bench]
        fn to_bytes(bencher: Bencher) {
            let a = sc_a();
            bencher.bench(|| black_box(&a).to_bytes());
        }

        #[divan::bench]
        fn from_bytes(bencher: Bencher) {
            let bytes = sc_a().to_bytes();
            bencher.bench(|| Scalar::from_bytes(black_box(&bytes)));
        }
    }

    /// Generate the group benchmarks, which are identical for G1 and G2.
    macro_rules! group_benches {
        ($name:ident, $group:path) => {
            mod $name {
                use super::*;
                use divan::{black_box, Bencher};
                use eccoxide::curve::group::CurveGroup;
                use $group::{Point, PointAffine};

                fn p() -> Point {
                    Point::GENERATOR
                }

                /// `2 * generator`, distinct from the generator.
                fn q() -> Point {
                    Point::GENERATOR.double()
                }

                #[divan::bench]
                fn add(bencher: Bencher) {
                    let (a, b) = (p(), q());
                    bencher.bench(|| black_box(&a) + black_box(&b));
                }

                #[divan::bench]
                fn double(bencher: Bencher) {
                    let a = p();
                    bencher.bench(|| black_box(&a).double());
                }

                /// Constant-time variable-base multiplication.
                #[divan::bench]
                fn scalar_mul(bencher: Bencher) {
                    let (a, k) = (p(), sc_a());
                    bencher.bench(|| black_box(&a) * black_box(&k));
                }

                /// Variable-time variable-base multiplication.
                #[divan::bench]
                fn scalar_mul_vartime(bencher: Bencher) {
                    let (a, k) = (p(), sc_a());
                    bencher.bench(|| black_box(&a).mul_vartime(black_box(&k)));
                }

                /// Fixed-base multiplication of the generator (comb table when
                /// the `table` feature is on).
                #[divan::bench]
                fn mul_base(bencher: Bencher) {
                    let k = sc_a();
                    bencher.bench(|| Point::mul_base(black_box(&k)));
                }

                #[divan::bench]
                fn to_affine(bencher: Bencher) {
                    // a genuine projective point (z != 1) so the inversion runs
                    let a = &p() + &q();
                    bencher.bench(|| black_box(&a).to_affine());
                }

                #[divan::bench]
                fn from_affine(bencher: Bencher) {
                    let a = PointAffine::GENERATOR;
                    bencher.bench(|| Point::from_affine(black_box(&a)));
                }

                #[divan::bench]
                fn compress(bencher: Bencher) {
                    let a = PointAffine::GENERATOR;
                    bencher.bench(|| black_box(&a).compress());
                }

                /// Prime-order-subgroup membership: the endomorphism check, to
                /// compare with the `scalar_mul` the definition would cost.
                #[divan::bench]
                fn is_in_subgroup(bencher: Bencher) {
                    let a = p();
                    bencher.bench(|| black_box(&a).is_in_subgroup());
                }

                /// Cofactor clearing: the map hash-to-curve ends with.
                #[divan::bench]
                fn clear_cofactor(bencher: Bencher) {
                    let a = p();
                    bencher.bench(|| black_box(&a).clear_cofactor());
                }

                /// RFC 9380 `hash_to_curve`: expand, two field elements, two
                /// SSWU maps plus isogenies, and the cofactor clearing.
                #[cfg(feature = "bls12-381-hash-to-curve")]
                #[divan::bench]
                fn hash_to_curve(bencher: Bencher) {
                    bencher.bench(|| {
                        Point::hash_to_curve(
                            black_box(b"benchmark message"),
                            black_box(b"QUUX-V01-CS02-with-eccoxide-bench"),
                        )
                    });
                }

                /// The non-uniform encoding: one field element and one map.
                #[cfg(feature = "bls12-381-hash-to-curve")]
                #[divan::bench]
                fn encode_to_curve(bencher: Bencher) {
                    bencher.bench(|| {
                        Point::encode_to_curve(
                            black_box(b"benchmark message"),
                            black_box(b"QUUX-V01-CS02-with-eccoxide-bench"),
                        )
                    });
                }

                #[divan::bench]
                fn decompress(bencher: Bencher) {
                    let a = PointAffine::GENERATOR;
                    let (x, sign) = a.compress();
                    let x = x.clone();
                    bencher.bench(|| PointAffine::decompress(black_box(&x), sign));
                }
            }
        };
    }

    group_benches!(g1, eccoxide::curve::bls12_381::g1);
    group_benches!(g2, eccoxide::curve::bls12_381::g2);

    // --- optimal-ate pairing ----------------------------------------------
    mod pairing {
        use super::*;
        use divan::{black_box, Bencher};
        use eccoxide::curve::bls12_381::pairing::{miller_loop, multi_miller_loop, pairing};
        use eccoxide::curve::bls12_381::{g1, g2};

        /// A non-generator G1 point, and a second one distinct from it.
        fn points(i: u64) -> (g1::PointAffine, g2::PointAffine) {
            let s = &sc_a() + &Scalar::from_u64(i);
            let t = &sc_b() + &Scalar::from_u64(i);
            (
                g1::Point::mul_base(&s).to_affine().unwrap(),
                g2::Point::mul_base(&t).to_affine().unwrap(),
            )
        }

        /// `e(P, Q)` for non-generator points: the Miller loop over `E'(Fp2)` in
        /// projective coordinates plus the easy/hard final exponentiation.
        #[divan::bench(sample_count = 20, sample_size = 1)]
        fn pairing_full(bencher: Bencher) {
            let (p, q) = points(0);
            bencher.bench(|| pairing(black_box(&p), black_box(&q)));
        }

        /// The Miller loop alone, without the final exponentiation.
        #[divan::bench(sample_count = 20, sample_size = 1)]
        fn miller(bencher: Bencher) {
            let (p, q) = points(0);
            bencher.bench(|| miller_loop(black_box(&p), black_box(&q)));
        }

        /// The final exponentiation alone.
        #[divan::bench(sample_count = 20, sample_size = 1)]
        fn final_exponentiation(bencher: Bencher) {
            let (p, q) = points(0);
            let m = miller_loop(&p, &q);
            bencher.bench(|| black_box(&m).final_exponentiation());
        }

        /// The shared Miller loop of a product of `N` pairings: compare with `N`
        /// times [`miller`] to see the saved `Fp12` squarings.
        #[divan::bench(sample_count = 20, sample_size = 1, consts = [2, 4])]
        fn multi_miller<const N: u64>(bencher: Bencher) {
            let pairs: Vec<_> = (0..N).map(points).collect();
            let terms: Vec<_> = pairs.iter().map(|(p, q)| (p, q)).collect();
            bencher.bench(|| multi_miller_loop(black_box(&terms)));
        }
    }
}
