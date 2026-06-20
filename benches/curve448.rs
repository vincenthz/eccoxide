//! Divan micro-benchmarks for curve448 ("Goldilocks").
//!
//! Covers the base field F(p), p = 2^448 - 2^224 - 1, and the x-only Montgomery
//! ladder (X448). There is no scalar field or Edwards group for this curve.
//!
//! Run with:
//!
//! ```text
//! cargo bench --bench curve448 --features curve448
//! cargo bench --bench curve448 --features curve448 -- field::mul
//! ```

fn main() {
    divan::main();
}

#[cfg(feature = "curve448")]
mod curve448 {
    use eccoxide::curve::curve448::*;

    // --- sample values ----------------------------------------------------

    /// A non-trivial, full-width field element.
    fn fe_a() -> FieldElement {
        let mut x = FieldElement::from_u64(0x0123_4567_89ab_cdef);
        for _ in 0..6 {
            x = x.square();
        }
        x
    }

    /// A second, distinct full-width field element.
    fn fe_b() -> FieldElement {
        let mut x = FieldElement::from_u64(0xfedc_ba98_7654_3210);
        for _ in 0..6 {
            x = x.square();
        }
        x
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
            bencher.bench(|| black_box(&a).to_bytes_le());
        }

        #[divan::bench]
        fn from_bytes(bencher: Bencher) {
            let bytes = fe_a().to_bytes_le();
            bencher.bench(|| FieldElement::from_bytes_le(black_box(&bytes)));
        }
    }

    // --- x-only Montgomery ladder -----------------------------------------
    mod montgomery {
        use super::*;
        use divan::{black_box, Bencher};

        // an arbitrary 56-byte (big-endian) scalar for the ladder
        fn scalar() -> [u8; 56] {
            core::array::from_fn(|i| (i as u8).wrapping_mul(7).wrapping_add(3))
        }

        #[divan::bench]
        fn scalar_mul(bencher: Bencher) {
            // the X448 ladder
            let (g, k) = (MontgomeryPoint::GENERATOR, scalar());
            bencher.bench(|| black_box(&g).scale_bytes(black_box(&k)));
        }

        #[divan::bench]
        fn u(bencher: Bencher) {
            let g = MontgomeryPoint::GENERATOR;
            bencher.bench(|| black_box(&g).u());
        }
    }
}
