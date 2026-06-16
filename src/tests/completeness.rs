//! Completeness checks for the projective `+` operator.
//!
//! The scalar-multiplication KATs only ever add *distinct* points, so on their
//! own they do not exercise the exceptional cases that the complete addition
//! formula is meant to handle. Since `add_or_double` no longer compares the two
//! operands and dispatches to a dedicated doubling routine (it relies on the
//! completeness of the addition formula instead), these cases are checked here
//! through the public `+` operator:
//!
//! * `P + P`  (must equal `2·P`)
//! * `P + (-P)` (must be the point at infinity)
//! * `P + O` and `O + P` (must be `P`)
//! * `(a·P) + (b·P) == (a + b)·P` (addition of two genuinely distinct points)

macro_rules! test_completeness {
    ($curve:ident) => {
        mod $curve {
            use crate::curve::sec2::$curve::{Point, Scalar};

            #[test]
            fn double_matches_scalar_mul() {
                // `P + P` goes through the complete addition formula (no
                // equality test / dedicated doubling), so it must agree with
                // `2·P` computed by scalar multiplication.
                let g = Point::generator();
                let sum = &g + &g;
                let twice = &g * &Scalar::from_u64(2);
                assert_eq!(sum, twice);
            }

            #[test]
            fn add_opposite_is_infinity() {
                let g = Point::generator();
                let inf = &g + &(-&g);
                assert_eq!(inf, Point::infinity());
            }

            #[test]
            fn add_infinity_is_identity() {
                let g = Point::generator();
                let inf = Point::infinity();
                assert_eq!(&g + &inf, g);
                assert_eq!(&inf + &g, g);
            }

            #[test]
            fn homomorphism() {
                // (a·P) + (b·P) == (a + b)·P, exercising the addition of two
                // genuinely distinct projective points.
                let g = Point::generator();
                let a = Scalar::from_u64(0x0123_4567_89ab_cdef);
                let b = Scalar::from_u64(0xfedc_ba98_7654_3210);
                let ab = &a + &b;
                let lhs = &(&g * &a) + &(&g * &b);
                let rhs = &g * &ab;
                assert_eq!(lhs, rhs);
            }

            #[test]
            fn ct_matches_vartime() {
                // The default `*` is the constant-time fixed-window routine; it
                // must agree with the variable-time double-and-add for every
                // scalar, including the edge cases (0, low/high window values).
                let g = Point::generator();
                // includes values that stress the wNAF recoding used by the
                // variable-time path: long runs of set bits (carry propagation),
                // the window boundaries and a top-bit-set value.
                for v in [
                    0u64,
                    1,
                    2,
                    15,
                    16,
                    17,
                    31,
                    32,
                    255,
                    256,
                    0x5555_5555_5555_5555,
                    0xaaaa_aaaa_aaaa_aaaa,
                    0xffff_ffff_ffff_ffff,
                    0x8000_0000_0000_0000,
                    0x0123_4567_89ab_cdef,
                ] {
                    let k = Scalar::from_u64(v);
                    assert_eq!(&g * &k, g.mul_vartime(&k), "mismatch for k = {}", v);
                }
                // a full-width scalar as well
                let mut big = Scalar::from_u64(0xfedc_ba98_7654_3210);
                for _ in 0..5 {
                    big = big.square();
                }
                assert_eq!(&g * &big, g.mul_vartime(&big));
            }

            #[test]
            fn mul_base_matches_generic() {
                // the fixed-base comb (mul_base) must agree with the general
                // variable-base scalar multiplication of the generator. This
                // validates the embedded precomputation tables end to end.
                let g = Point::generator();
                for v in [0u64, 1, 2, 15, 16, 17, 255, 256, 0x0123_4567_89ab_cdef] {
                    let k = Scalar::from_u64(v);
                    assert_eq!(Point::mul_base(&k), &g * &k, "mul_base mismatch for k = {}", v);
                }
                // a full-width scalar as well
                let mut big = Scalar::from_u64(0xfedc_ba98_7654_3210);
                for _ in 0..5 {
                    big = big.square();
                }
                assert_eq!(Point::mul_base(&big), &g * &big);
            }
        }
    };
}

#[cfg(feature = "p192k1")]
test_completeness!(p192k1);
#[cfg(feature = "p192r1")]
test_completeness!(p192r1);
#[cfg(feature = "p224k1")]
test_completeness!(p224k1);
#[cfg(feature = "p224r1")]
test_completeness!(p224r1);
#[cfg(feature = "p256k1")]
test_completeness!(p256k1);
#[cfg(feature = "p256r1")]
test_completeness!(p256r1);
#[cfg(feature = "p384r1")]
test_completeness!(p384r1);
#[cfg(feature = "p521r1")]
test_completeness!(p521r1);
