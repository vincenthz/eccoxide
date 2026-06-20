//! Divan micro-benchmarks for curve25519 / edwards25519.
//!
//! Covers the base field F(p), the scalar field F(l), the twisted-Edwards group
//! and the x-only Montgomery ladder.
//!
//! Run with:
//!
//! ```text
//! cargo bench --bench curve25519
//! cargo bench --bench curve25519 -- edwards     # only the Edwards group
//! cargo bench --bench curve25519 -- field::mul
//! ```

fn main() {
    divan::main();
}

#[cfg(feature = "curve25519")]
mod curve25519 {
    use eccoxide::curve::curve25519::*;

    // --- sample values ----------------------------------------------------

    /// A non-trivial, full-width field element (the generator's x).
    fn fe_a() -> FieldElement {
        Point::GENERATOR.to_affine().0
    }

    /// A second, distinct full-width field element (the generator's y).
    fn fe_b() -> FieldElement {
        Point::GENERATOR.to_affine().1
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

    // --- prime field F(p) -------------------------------------------------
    mod field {
        use super::*;
        use divan::{black_box, Bencher};

        #[divan::bench]
        fn add(bencher: Bencher) {
            let (a, b) = (fe_a(), fe_b());
            bencher.bench(|| black_box(&a) + black_box(&b));
        }

        #[divan::bench]
        fn sub(bencher: Bencher) {
            let (a, b) = (fe_a(), fe_b());
            bencher.bench(|| black_box(&a) - black_box(&b));
        }

        #[divan::bench]
        fn mul(bencher: Bencher) {
            let (a, b) = (fe_a(), fe_b());
            bencher.bench(|| black_box(&a) * black_box(&b));
        }

        #[divan::bench]
        fn square(bencher: Bencher) {
            let a = fe_a();
            bencher.bench(|| black_box(&a).square());
        }

        #[divan::bench]
        fn double(bencher: Bencher) {
            let a = fe_a();
            bencher.bench(|| black_box(&a).double());
        }

        #[divan::bench]
        fn neg(bencher: Bencher) {
            let a = fe_a();
            bencher.bench(|| -black_box(&a));
        }

        #[divan::bench]
        fn inverse(bencher: Bencher) {
            let a = fe_a();
            bencher.bench(|| black_box(&a).inverse());
        }

        #[divan::bench]
        fn sqrt(bencher: Bencher) {
            // square first so the input is guaranteed to be a residue
            let a = fe_a().square();
            bencher.bench(|| black_box(&a).sqrt());
        }

        #[divan::bench]
        fn to_bytes(bencher: Bencher) {
            let a = fe_a();
            bencher.bench(|| black_box(&a).to_bytes());
        }

        #[divan::bench]
        fn from_bytes(bencher: Bencher) {
            let bytes = fe_a().to_bytes();
            bencher.bench(|| FieldElement::from_bytes(black_box(&bytes)));
        }
    }

    // --- scalar field F(l) ------------------------------------------------
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
        fn inverse_safegcd(bencher: Bencher) {
            let a = sc_a();
            bencher.bench(|| black_box(&a).inverse_safegcd());
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

    // --- twisted-Edwards group --------------------------------------------
    mod edwards {
        use super::*;
        use divan::{black_box, Bencher};

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

        #[divan::bench]
        fn scalar_mul(bencher: Bencher) {
            // constant-time double-and-add over the complete addition formula
            let (a, k) = (p(), sc_a());
            bencher.bench(|| black_box(&a) * black_box(&k));
        }

        #[divan::bench]
        fn to_affine(bencher: Bencher) {
            // a genuine projective point (z != 1) so the modular inversion runs
            let a = &p() + &q();
            bencher.bench(|| black_box(&a).to_affine());
        }

        #[divan::bench]
        fn compress(bencher: Bencher) {
            let a = p();
            bencher.bench(|| black_box(&a).compress());
        }

        #[divan::bench]
        fn decompress(bencher: Bencher) {
            let (y, sign) = p().compress();
            bencher.bench(|| Point::decompress(black_box(&y), sign));
        }

        #[divan::bench]
        fn to_montgomery(bencher: Bencher) {
            let a = p();
            bencher.bench(|| black_box(&a).to_montgomery());
        }
    }

    // --- x-only Montgomery ladder -----------------------------------------
    mod montgomery {
        use super::*;
        use divan::{black_box, Bencher};

        #[divan::bench]
        fn scalar_mul(bencher: Bencher) {
            // the X25519 ladder
            let (g, k) = (MontgomeryPoint::GENERATOR, sc_a());
            bencher.bench(|| black_box(&g).scale(black_box(&k)));
        }

        #[divan::bench]
        fn u(bencher: Bencher) {
            let g = MontgomeryPoint::GENERATOR;
            bencher.bench(|| black_box(&g).u());
        }
    }
}
