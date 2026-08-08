//! BLS12-381 G1
//!
//! G1 is the prime-order-`r` subgroup of the short Weierstrass curve
//! `E(Fp): y^2 = x^3 + 4` (so `a = 0`, `b = 4`) defined over the base field
//! [`Fp`](super::fp::Fp). Its scalar field is [`Scalar`](super::scalar::Scalar).
//!
//! This wires the generic curve machinery (affine / projective points, the
//! complete addition formulas and the fixed-base comb) into the [`CurveGroup`]
//! abstraction.
//!
//! # G1 versus the curve
//!
//! `E(Fp)` has order `h·r`, with the 126-bit cofactor `h = (x - 1)²/3` for the
//! seed `x`, so satisfying the curve equation is *weaker* than being in the
//! prime-order group G1 the pairing and the protocols above it are defined on.
//! [`Point::is_in_subgroup`] tests the difference, and
//! [`Point::clear_cofactor`] establishes it.
//!
//! [`CurveGroup`]: crate::curve::group::CurveGroup

use crate::curve::bls12_381::fp::Fp as FieldElement;
#[cfg(feature = "bls12-381-hash-to-curve")]
use crate::curve::bls12_381::hash_to_curve;
use crate::curve::bls12_381::scalar::Scalar;
use crate::curve::bls12_381::BLS_X;
use crate::curve::field::{Field, Sign};
use crate::curve::{
    affine, projective,
    weierstrass::{WeierstrassCurve, WeierstrassCurveA0},
};
use crate::mp::ct::{Choice, CtEqual};
use crate::params::bls12_381::*;
use crate::{
    bls12_381_define_point_serialization, fiat_define_weierstrass_curve,
    fiat_define_weierstrass_curve_a0, fiat_define_weierstrass_points,
};

fiat_define_weierstrass_curve!(FieldElement);
fiat_define_weierstrass_curve_a0!(FieldElement);
fiat_define_weierstrass_points!(FieldElement);
bls12_381_define_point_serialization!(FieldElement);

/// `β`, a primitive cube root of unity in `Fp`; see [`BETA_BYTES`].
const BETA: FieldElement = FieldElement::from_bytes_unchecked(&BETA_BYTES);

impl Point {
    /// The endomorphism `σ(x, y) = (β·x, y)` of the curve, where `β` is a
    /// primitive cube root of unity in `Fp`.
    ///
    /// It is defined on all of `E(Fp)` — `(β·x)³ + 4 = x³ + 4` since `β³ = 1` —
    /// and on G1 it is multiplication by the scalar `-x²`, one of the two
    /// primitive cube roots of unity modulo `r`. That is what makes it both a
    /// GLV endomorphism for scalar multiplication and, in the other direction,
    /// the cheap witness [`Self::is_in_subgroup`] uses.
    fn sigma(&self) -> Self {
        Point(projective::Point {
            x: &self.0.x * &BETA,
            y: self.0.y.clone(),
            z: self.0.z.clone(),
        })
    }

    /// `[|x|]·self`, for the magnitude of the seed `x`, by double-and-add.
    ///
    /// `|x|` is a public 64-bit constant with only 6 bits set, so this is 63
    /// doublings and 5 additions, and following its bits leaks nothing about
    /// the point. The additions use the complete formula, so no case is special.
    fn mul_by_abs_x(&self) -> Self {
        // the top bit of |x| is set, and consumed by starting the loop from self
        let mut acc = self.0.clone();
        for i in (0..63).rev() {
            acc = acc.double_a0::<Curve>();
            if (BLS_X >> i) & 1 == 1 {
                acc = acc.add_different_a0::<Curve>(&self.0);
            }
        }
        Point(acc)
    }

    /// Whether the point lies in the prime-order-`r` subgroup G1, and not merely
    /// on the curve.
    ///
    /// Points read from an untrusted source have to be checked: `E(Fp)` has the
    /// 126-bit cofactor `h = (x - 1)²/3`, and a point of order `h·r` (or of any
    /// order sharing a factor with `h`) breaks the assumptions the pairing-based
    /// protocols rest on, even though it satisfies the curve equation. The
    /// standard decoders ([`PointAffine::from_compressed`] and friends) already
    /// run this check.
    ///
    /// The test is the endomorphism one of section 4 of
    /// <https://eprint.iacr.org/2021/1130> (proof of correctness revised in
    /// <https://eprint.iacr.org/2022/352>): for `P` on the curve,
    ///
    /// ```text
    /// P ∈ G1   ⟺   σ(P) = [-x²]P
    /// ```
    ///
    /// where `σ(x, y) = (β·x, y)` for a primitive cube root of unity `β`. Since
    /// the seed is only 64 bits, evaluating the right-hand side costs two
    /// multiplications by it and nothing else, about 2.5x less than the
    /// multiplication by the 255-bit `r` the definition `[r]P == 𝒪` would ask
    /// for.
    ///
    /// The identity is in the subgroup, so this is true of it.
    pub fn is_in_subgroup(&self) -> Choice {
        // x² = |x|², the sign of the seed being irrelevant to a square
        let minus_x2_p = -&self.mul_by_abs_x().mul_by_abs_x();
        self.sigma().0.ct_eq(&minus_x2_p.0)
    }

    /// Send an arbitrary point of the curve into the prime-order subgroup G1,
    /// by multiplying it by the effective cofactor `1 - x = 0xd201000000010001`.
    ///
    /// This is the last step of hash-to-curve: a map onto `E(Fp)` lands anywhere
    /// on the curve, and this brings the result into the group the protocols
    /// work in. The output is in G1 for *any* input, the identity included, and
    /// the map is surjective onto G1.
    ///
    /// Multiplying by `1 - x` rather than by the 126-bit cofactor `h = (x-1)²/3`
    /// is enough because `E(Fp)` is not cyclic — its `h`-torsion splits so that
    /// the short factor already kills it (section 5 of
    /// <https://eprint.iacr.org/2019/403>). `1 - x` is exactly the `h_eff` that
    /// section 8.8.1 of RFC 9380 prescribes.
    ///
    /// Note this is not the identity on points that are already in G1: there it
    /// is multiplication by the scalar `1 - x`. Use [`Self::is_in_subgroup`] to
    /// *test* membership; this one *establishes* it.
    ///
    /// The result is projective because clearing the cofactor of a point of
    /// order dividing `h` gives the identity, which the affine type cannot hold.
    pub fn clear_cofactor(&self) -> Self {
        // [1 - x]P = P + [|x|]P, the seed being negative
        self + &self.mul_by_abs_x()
    }
}

impl PointAffine {
    /// Whether the point lies in the prime-order-`r` subgroup G1, and not merely
    /// on the curve; see [`Point::is_in_subgroup`], which this defers to.
    pub fn is_in_subgroup(&self) -> Choice {
        Point::from(self).is_in_subgroup()
    }

    /// Send an arbitrary point of the curve into the prime-order subgroup G1;
    /// see [`Point::clear_cofactor`], which this defers to.
    ///
    /// The result is projective: the cleared point can be the identity.
    pub fn clear_cofactor(&self) -> Point {
        Point::from(self).clear_cofactor()
    }
}

#[cfg(feature = "bls12-381-hash-to-curve")]
impl Point {
    /// Wrap a bare projective point, for the hash-to-curve tests: the maps
    /// build their result coordinate by coordinate, from outside this module.
    #[cfg(test)]
    pub(crate) fn from_projective(p: projective::Point<FieldElement>) -> Self {
        Point(p)
    }

    /// Hash a message to a point of G1, as the RFC 9380 suite
    /// `BLS12381G1_XMD:SHA-256_SSWU_RO_`.
    ///
    /// The result is uniformly distributed over G1 and the map is
    /// indifferentiable from a random oracle, which is what protocols building
    /// on it (BLS signatures, verifiable random functions, commitments to
    /// arbitrary data) assume. Prefer this over
    /// [`Self::encode_to_curve`] unless the input is already uniform.
    ///
    /// `dst` is the domain separation tag: pick one string per application and
    /// per use within it, and keep it fixed. RFC 9380 section 3.1 recommends
    /// embedding the suite name, e.g.
    /// `b"MYPROTOCOL-V01-CS01-with-BLS12381G1_XMD:SHA-256_SSWU_RO_"`. A tag longer
    /// than 255 bytes is hashed down as the specification prescribes, so any
    /// length works.
    ///
    /// The message is hashed to two field elements, each mapped to the curve
    /// and added, and the sum has its cofactor cleared — see the
    /// [`hash_to_curve` module](super::hash_to_curve) for the pipeline.
    pub fn hash_to_curve(msg: &[u8], dst: &[u8]) -> Self {
        let [u0, u1] = hash_to_curve::hash_to_field_g1(msg, dst);
        let q0 = Point(hash_to_curve::map_to_curve_g1(&u0));
        let q1 = Point(hash_to_curve::map_to_curve_g1(&u1));
        (&q0 + &q1).clear_cofactor()
    }

    /// Encode a message to a point of G1, as the RFC 9380 suite
    /// `BLS12381G1_XMD:SHA-256_SSWU_NU_`.
    ///
    /// Cheaper than [`Self::hash_to_curve`] — one field element and one map
    /// instead of two — but *non-uniform*: its image is a fraction of G1,
    /// and it is not indifferentiable from a random oracle. Only use it when
    /// the caller already guarantees a uniform input, or when the protocol
    /// explicitly asks for the `NU_` suite.
    ///
    /// See [`Self::hash_to_curve`] about the domain separation tag.
    pub fn encode_to_curve(msg: &[u8], dst: &[u8]) -> Self {
        let u = hash_to_curve::encode_to_field_g1(msg, dst);
        Point(hash_to_curve::map_to_curve_g1(&u)).clear_cofactor()
    }
}

#[cfg(test)]
mod tests {
    use super::{FieldElement, Point, PointAffine, Scalar};
    use crate::curve::group::CurveGroup;

    fn s(v: u64) -> Scalar {
        Scalar::from_u64(v)
    }

    #[test]
    fn generator_on_curve() {
        // the generator coordinates must satisfy the curve equation
        let (gx, gy) = PointAffine::GENERATOR.to_coordinate();
        assert!(PointAffine::from_coordinate(gx, gy).is_some());
    }

    #[test]
    fn identity_laws() {
        let g = Point::GENERATOR;
        let id = Point::INFINITY;
        assert_eq!(&g + &id, g);
        assert_eq!(&id + &g, g);
        assert_eq!(&g - &g, id);
        assert_eq!(&g + (-&g), id);
    }

    #[test]
    fn double_is_add() {
        let g = Point::GENERATOR;
        assert_eq!(g.double(), &g + &g);
    }

    #[test]
    fn scalar_mul_matches_addition() {
        let g = Point::GENERATOR;
        let g2 = g.double();
        let g3 = &g2 + &g;
        let three = s(3);
        assert_eq!(&g * &three, g3);
        assert_eq!(Point::mul_base(&three), g3);
        assert_eq!(g.mul_vartime(&three), g3);
    }

    #[test]
    fn mul_base_matches_variable_base() {
        // the fixed-base comb table must agree with variable-base multiplication
        let g = Point::GENERATOR;
        for v in [1u64, 2, 7, 255, 256, 0x1_0001, 0x1234_5678_9abc_def0] {
            let k = s(v);
            assert_eq!(Point::mul_base(&k), &g * &k, "mul_base != g*k for {}", v);
            assert_eq!(g.mul_vartime(&k), &g * &k, "mul_vartime != g*k for {}", v);
        }
    }

    #[test]
    fn scalar_mul_is_linear() {
        let g = Point::GENERATOR;
        let a = Scalar::from(0x1234_5678_9abc_def0u64);
        let b = Scalar::from(0x0fed_cba9_8765_4321u64);
        let ab = &a + &b;
        assert_eq!(
            Point::mul_base(&ab),
            Point::mul_base(&a) + Point::mul_base(&b)
        );
        assert_eq!(&g * &ab, &g * &a + &g * &b);
    }

    #[test]
    fn compress_decompress_roundtrip() {
        let p = (Point::GENERATOR * &s(12345)).to_affine().unwrap();
        let (x, sign) = p.compress();
        let x = x.clone();
        let recovered = PointAffine::decompress(&x, sign).into_option().unwrap();
        assert_eq!(recovered, p);
    }

    #[test]
    fn to_affine_from_affine_roundtrip() {
        let p = Point::GENERATOR * &s(98765);
        let affine = p.to_affine().unwrap();
        assert_eq!(Point::from_affine(&affine), p);
    }

    #[test]
    fn generator_field_element_bytes_roundtrip() {
        let (gx, _) = PointAffine::GENERATOR.to_coordinate();
        let bytes = gx.to_bytes();
        assert_eq!(&FieldElement::from_bytes(&bytes).unwrap(), gx);
    }

    crate::bls12_381_point_serialization_unittest!();

    /// Prime-order-subgroup membership, and the endomorphism it is built on.
    mod subgroup {
        use super::super::{Curve, Point, PointAffine};
        use super::{s, FieldElement, Scalar};
        use crate::curve::bls12_381::BLS_X;
        use crate::curve::field::Sign;
        use crate::params::bls12_381::ORDER_BYTES;

        /// Points of the curve that are **not** in the prime-order subgroup, as their
        /// `(compressed, uncompressed)` encodings.
        ///
        /// The first two are the decompressions of the smallest x-coordinates that are
        /// on the curve at all, `x = 4` and `x = 5`, with opposite y-signs so that the
        /// sort flag is pinned both ways; the third is `G` plus the first, the shape an
        /// attacker sends — a genuine G1 point with a torsion component added.
        ///
        /// Generated, and checked against `[r]P != 𝒪`, by `sage/bls12_381.sage`.
        const OFF_SUBGROUP: &[([u8; 48], [u8; 96])] = &[
            (
                [
                    0x80, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04,
                ],
                [
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04, 0x0a, 0x98, 0x9b, 0xad,
                    0xd4, 0x0d, 0x62, 0x12, 0xb3, 0x3c, 0xff, 0xc3, 0xf3, 0x76, 0x3e, 0x9b, 0xc7,
                    0x60, 0xf9, 0x88, 0xc9, 0x92, 0x6b, 0x26, 0xda, 0x9d, 0xd8, 0x5e, 0x92, 0x84,
                    0x83, 0x44, 0x63, 0x46, 0xb8, 0xed, 0x00, 0xe1, 0xde, 0x5d, 0x5e, 0xa9, 0x3e,
                    0x35, 0x4a, 0xbe, 0x70, 0x6c,
                ],
            ),
            (
                [
                    0xa0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x05,
                ],
                [
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x05, 0x0d, 0x3c, 0x6d, 0xa1,
                    0x21, 0x1e, 0xbe, 0x79, 0x7b, 0xc0, 0x79, 0x0f, 0x1e, 0x6e, 0x7d, 0x66, 0x9b,
                    0x18, 0x0a, 0x8e, 0x59, 0x19, 0x68, 0x25, 0x50, 0x6d, 0x2b, 0xb2, 0x18, 0x5f,
                    0x53, 0x71, 0x5d, 0xf0, 0x92, 0xc8, 0xa7, 0xce, 0xb6, 0x48, 0x43, 0xea, 0x7d,
                    0xf6, 0x7d, 0xba, 0xd6, 0x0d,
                ],
            ),
            (
                [
                    0x97, 0xbc, 0xbb, 0xfd, 0xd2, 0x44, 0x2c, 0x32, 0x81, 0x50, 0xf6, 0x54, 0x65,
                    0xbd, 0x7b, 0x9c, 0x4f, 0xf3, 0x6e, 0x35, 0x26, 0x1a, 0xd3, 0x54, 0x92, 0x22,
                    0xe5, 0x32, 0x75, 0x8a, 0x1c, 0xf0, 0x94, 0x5b, 0xa1, 0x33, 0xec, 0x51, 0x35,
                    0x17, 0xb4, 0xea, 0x9d, 0xe0, 0x98, 0xa0, 0x37, 0xf9,
                ],
                [
                    0x17, 0xbc, 0xbb, 0xfd, 0xd2, 0x44, 0x2c, 0x32, 0x81, 0x50, 0xf6, 0x54, 0x65,
                    0xbd, 0x7b, 0x9c, 0x4f, 0xf3, 0x6e, 0x35, 0x26, 0x1a, 0xd3, 0x54, 0x92, 0x22,
                    0xe5, 0x32, 0x75, 0x8a, 0x1c, 0xf0, 0x94, 0x5b, 0xa1, 0x33, 0xec, 0x51, 0x35,
                    0x17, 0xb4, 0xea, 0x9d, 0xe0, 0x98, 0xa0, 0x37, 0xf9, 0x06, 0xd1, 0xd4, 0xf6,
                    0x58, 0x0f, 0x49, 0xb4, 0xe0, 0xa9, 0x85, 0x09, 0xff, 0xd1, 0x8f, 0x24, 0xaf,
                    0xca, 0xda, 0x36, 0xfd, 0x0d, 0x44, 0xe9, 0xfc, 0x9e, 0x5f, 0x0c, 0x19, 0xdf,
                    0x3e, 0xc0, 0x14, 0x74, 0xee, 0xfc, 0x65, 0x9d, 0x57, 0xd1, 0x49, 0xb9, 0x7c,
                    0xa8, 0x99, 0x01, 0x0a, 0x5d,
                ],
            ),
        ];

        /// `[r]p`, the definition the fast check has to agree with: a point of
        /// the curve is in the prime-order subgroup exactly when the group order
        /// kills it. Deliberately computed the slow way, by a generic scalar
        /// multiplication by `r` — which no `Scalar` could hold, `r` being the
        /// modulus of that very field.
        fn mul_by_order(p: &Point) -> Point {
            Point(p.0.scale_a0::<Curve>(&ORDER_BYTES))
        }

        /// Points of the curve found by decompressing small x-coordinates. The
        /// cofactor of `E(Fp)` is 126 bits wide, so none of these is in G1 by
        /// any plausible chance, and about half the x tried are on the curve.
        fn small_x_points() -> Vec<PointAffine> {
            let points: Vec<_> = (1..32u64)
                .filter_map(|i| {
                    PointAffine::decompress(&FieldElement::from(i), Sign::Positive).into_option()
                })
                .collect();
            assert!(!points.is_empty(), "no small x gave a point of the curve");
            points
        }

        /// Multiples of the generator are in G1, whichever representation they
        /// are held in, and so is the identity.
        #[test]
        fn generator_multiples_are_in_the_subgroup() {
            for v in [1u64, 2, 3, 12345, 0xdead_beef] {
                let p = Point::GENERATOR * &s(v);
                assert!(p.is_in_subgroup().is_true(), "[{}]G rejected", v);
                assert!(p.to_affine().unwrap().is_in_subgroup().is_true());
            }
            assert!(Point::INFINITY.is_in_subgroup().is_true());
        }

        /// The check must agree with `[r]P == 𝒪` on points of the curve, both in
        /// the subgroup and outside it — that equivalence is the whole claim of
        /// the endomorphism test.
        #[test]
        fn matches_multiplication_by_the_order() {
            let mut points: Vec<Point> = [1u64, 3, 12345]
                .iter()
                .map(|v| Point::GENERATOR * &s(*v))
                .collect();
            points.push(Point::INFINITY);
            points.extend(small_x_points().iter().map(Point::from));
            // the shape an attacker would send: a point of G1 with a torsion
            // component added, so the order is a proper multiple of r
            points.extend(
                small_x_points()
                    .iter()
                    .map(|q| &Point::GENERATOR + &Point::from(q)),
            );

            let (mut inside, mut outside) = (0, 0);
            for p in points {
                let expected = mul_by_order(&p) == Point::INFINITY;
                assert_eq!(
                    p.is_in_subgroup().is_true(),
                    expected,
                    "the endomorphism check disagrees with [r]P == O"
                );
                if expected {
                    inside += 1;
                } else {
                    outside += 1;
                }
            }
            assert!(inside > 0 && outside > 0, "one side was never exercised");
        }

        /// The eigenvalue of `σ` on G1, which is what pins `β`: the other
        /// non-trivial cube root of unity would give `[x² - 1]` instead.
        #[test]
        fn sigma_is_multiplication_by_minus_x_squared() {
            let abs_x = Scalar::from_u64(BLS_X);
            let lambda = -&(&abs_x * &abs_x);
            for v in [1u64, 2, 12345] {
                let p = Point::GENERATOR * &s(v);
                assert_eq!(p.sigma(), &p * &lambda, "sigma != [-x^2] on [{}]G", v);
            }
        }

        /// `σ` is an endomorphism of the whole curve, not just of G1.
        #[test]
        fn sigma_stays_on_the_curve() {
            for p in small_x_points() {
                let q = Point::from(&p).sigma().to_affine().unwrap();
                let (x, y) = q.to_coordinate();
                assert!(PointAffine::from_coordinate(x, y).is_some());
            }
        }

        /// The double-and-add over the seed agrees with a generic scalar
        /// multiplication by it.
        #[test]
        fn mul_by_abs_x_matches_scalar_mul() {
            let k = Scalar::from_u64(BLS_X);
            for v in [1u64, 7, 12345] {
                let p = Point::GENERATOR * &s(v);
                assert_eq!(p.mul_by_abs_x(), &p * &k);
            }
        }

        /// The hardcoded vectors are what they claim to be — on the curve, and
        /// outside the subgroup by the definition — so the standard decoders
        /// must refuse them and the `_oncurve_only` ones must accept them.
        ///
        /// Unlike [`small_x_points`] these are fixed bytes, so they also pin the
        /// encodings themselves: a point that stops being rejected, or an
        /// encoding that stops round-tripping, shows up here.
        #[test]
        fn off_subgroup_kat() {
            for (i, (compressed, uncompressed)) in OFF_SUBGROUP.iter().enumerate() {
                let p = PointAffine::from_compressed_oncurve_only(compressed)
                    .unwrap_or_else(|| panic!("vector {} is not a point of the curve", i));
                let q = PointAffine::from_uncompressed_oncurve_only(uncompressed)
                    .unwrap_or_else(|| panic!("vector {} is not a point of the curve", i));
                assert_eq!(p, q, "the two encodings of vector {} disagree", i);

                // outside the subgroup by the definition, and by the fast check
                assert_ne!(
                    mul_by_order(&Point::from(&p)),
                    Point::INFINITY,
                    "vector {} is in the subgroup after all",
                    i
                );
                assert!(!p.is_in_subgroup().is_true(), "vector {} accepted", i);
                assert!(
                    !Point::from(&p).is_in_subgroup().is_true(),
                    "vector {} accepted",
                    i
                );

                // so the checking decoders reject both encodings, for both types
                assert!(PointAffine::from_compressed(compressed).is_none());
                assert!(Point::from_compressed(compressed).is_none());
                assert!(PointAffine::from_uncompressed(uncompressed).is_none());
                assert!(Point::from_uncompressed(uncompressed).is_none());

                // and the vectors are the canonical encodings of that point
                assert_eq!(
                    &p.to_compressed(),
                    compressed,
                    "vector {} re-encodes differently",
                    i
                );
                assert_eq!(
                    &p.to_uncompressed(),
                    uncompressed,
                    "vector {} re-encodes differently",
                    i
                );
            }
        }

        /// Cofactor clearing lands in the subgroup whatever the input point is:
        /// multiples of the generator, points off the subgroup, sums of the two,
        /// the identity, and the hardcoded vectors.
        #[test]
        fn clear_cofactor_lands_in_the_subgroup() {
            let mut points: Vec<Point> = [1u64, 3, 12345]
                .iter()
                .map(|v| Point::GENERATOR * &s(*v))
                .collect();
            points.push(Point::INFINITY);
            points.extend(small_x_points().iter().map(Point::from));
            points.extend(
                small_x_points()
                    .iter()
                    .map(|q| &Point::GENERATOR + &Point::from(q)),
            );
            points.extend(
                OFF_SUBGROUP.iter().map(|(c, _)| {
                    Point::from(&PointAffine::from_compressed_oncurve_only(c).unwrap())
                }),
            );

            for p in points {
                let cleared = p.clear_cofactor();
                assert!(cleared.is_in_subgroup().is_true(), "cleared point rejected");
                assert_eq!(
                    mul_by_order(&cleared),
                    Point::INFINITY,
                    "[r] does not kill the cleared point"
                );
            }
        }

        /// It is multiplication by the `h_eff` of RFC 9380 section 8.8.1,
        /// `1 - x`, and not by anything else — checked against a generic
        /// multiplication by those bytes, for points on and off the subgroup.
        #[test]
        fn clear_cofactor_is_multiplication_by_h_eff() {
            // 1 - x = 1 + |x|, the seed being negative
            const H_EFF: u64 = BLS_X + 1;
            assert_eq!(H_EFF, 0xd201_0000_0001_0001);

            let mut points: Vec<Point> = vec![Point::GENERATOR, Point::GENERATOR * &s(12345)];
            points.extend(small_x_points().iter().map(Point::from));

            for p in points {
                let expected = Point(p.0.scale_a0::<Curve>(&H_EFF.to_be_bytes()));
                assert_eq!(p.clear_cofactor(), expected);
            }
            // on the subgroup that is just the scalar 1 - x, so it is *not* a no-op
            let g = Point::GENERATOR;
            assert_eq!(g.clear_cofactor(), &g * &Scalar::from_u64(H_EFF));
            assert_ne!(g.clear_cofactor(), g);
        }

        /// The identity clears to the identity, a stray point does not collapse
        /// to it, and the affine entry point agrees with the projective one.
        #[test]
        fn clear_cofactor_edge_cases() {
            assert_eq!(Point::INFINITY.clear_cofactor(), Point::INFINITY);
            for p in small_x_points() {
                assert_ne!(
                    Point::from(&p).clear_cofactor(),
                    Point::INFINITY,
                    "a curve point cleared to the identity"
                );
                assert_eq!(p.clear_cofactor(), Point::from(&p).clear_cofactor());
            }
        }
    }

    /// Known-answer vectors for the standard G1 encodings: `k * G` for small
    /// `k`, byte for byte.
    ///
    /// `k = 1` is the published encoding of the generator: its x-coordinate
    /// `0x17f1d3a7...` with the compression flag set, and the sort flag clear
    /// since its y is below `(p - 1)/2`. `k = 2` and `k = 0xdeadbeef` have the
    /// sort flag set instead, so the vectors pin the flag both ways. The
    /// multiples were computed independently.
    #[test]
    fn serialization_kat() {
        const COMPRESSED: &[(u64, [u8; 48])] = &[
            (
                1,
                [
                    0x97, 0xf1, 0xd3, 0xa7, 0x31, 0x97, 0xd7, 0x94, 0x26, 0x95, 0x63, 0x8c, 0x4f,
                    0xa9, 0xac, 0x0f, 0xc3, 0x68, 0x8c, 0x4f, 0x97, 0x74, 0xb9, 0x05, 0xa1, 0x4e,
                    0x3a, 0x3f, 0x17, 0x1b, 0xac, 0x58, 0x6c, 0x55, 0xe8, 0x3f, 0xf9, 0x7a, 0x1a,
                    0xef, 0xfb, 0x3a, 0xf0, 0x0a, 0xdb, 0x22, 0xc6, 0xbb,
                ],
            ),
            (
                2,
                [
                    0xa5, 0x72, 0xcb, 0xea, 0x90, 0x4d, 0x67, 0x46, 0x88, 0x08, 0xc8, 0xeb, 0x50,
                    0xa9, 0x45, 0x0c, 0x97, 0x21, 0xdb, 0x30, 0x91, 0x28, 0x01, 0x25, 0x43, 0x90,
                    0x2d, 0x0a, 0xc3, 0x58, 0xa6, 0x2a, 0xe2, 0x8f, 0x75, 0xbb, 0x8f, 0x1c, 0x7c,
                    0x42, 0xc3, 0x9a, 0x8c, 0x55, 0x29, 0xbf, 0x0f, 0x4e,
                ],
            ),
            (
                3,
                [
                    0x89, 0xec, 0xe3, 0x08, 0xf9, 0xd1, 0xf0, 0x13, 0x17, 0x65, 0x21, 0x2d, 0xec,
                    0xa9, 0x96, 0x97, 0xb1, 0x12, 0xd6, 0x1f, 0x9b, 0xe9, 0xa5, 0xf1, 0xf3, 0x78,
                    0x0a, 0x51, 0x33, 0x5b, 0x3f, 0xf9, 0x81, 0x74, 0x7a, 0x0b, 0x2c, 0xa2, 0x17,
                    0x9b, 0x96, 0xd2, 0xc0, 0xc9, 0x02, 0x4e, 0x52, 0x24,
                ],
            ),
            (
                12345,
                [
                    0x85, 0x30, 0xc1, 0xbd, 0xc4, 0xcd, 0x6b, 0x14, 0x08, 0xbe, 0x09, 0x33, 0xc4,
                    0xa4, 0x1a, 0xc3, 0x51, 0x33, 0x50, 0xee, 0xf3, 0x68, 0x50, 0xb8, 0x04, 0x70,
                    0x8e, 0x1f, 0x33, 0x89, 0x32, 0xce, 0x01, 0xb6, 0x55, 0xa1, 0x63, 0x34, 0x4a,
                    0x45, 0x00, 0xb2, 0x81, 0xc8, 0x75, 0x0c, 0x46, 0x1f,
                ],
            ),
            (
                0xdead_beef,
                [
                    0xac, 0xcc, 0xb5, 0xba, 0xb2, 0x94, 0x4a, 0x1b, 0xdc, 0x72, 0x1c, 0x97, 0xf3,
                    0xaf, 0xfa, 0x03, 0x5d, 0x50, 0x7c, 0x78, 0xfe, 0x44, 0x2a, 0x92, 0x84, 0x98,
                    0x2b, 0xd4, 0xc2, 0x76, 0x17, 0xb3, 0x3f, 0x1d, 0x46, 0xe8, 0x19, 0x1a, 0x1e,
                    0xda, 0x03, 0xd7, 0x3c, 0x35, 0x77, 0x52, 0xd2, 0x19,
                ],
            ),
        ];
        const UNCOMPRESSED: &[(u64, [u8; 96])] = &[
            (
                1,
                [
                    0x17, 0xf1, 0xd3, 0xa7, 0x31, 0x97, 0xd7, 0x94, 0x26, 0x95, 0x63, 0x8c, 0x4f,
                    0xa9, 0xac, 0x0f, 0xc3, 0x68, 0x8c, 0x4f, 0x97, 0x74, 0xb9, 0x05, 0xa1, 0x4e,
                    0x3a, 0x3f, 0x17, 0x1b, 0xac, 0x58, 0x6c, 0x55, 0xe8, 0x3f, 0xf9, 0x7a, 0x1a,
                    0xef, 0xfb, 0x3a, 0xf0, 0x0a, 0xdb, 0x22, 0xc6, 0xbb, 0x08, 0xb3, 0xf4, 0x81,
                    0xe3, 0xaa, 0xa0, 0xf1, 0xa0, 0x9e, 0x30, 0xed, 0x74, 0x1d, 0x8a, 0xe4, 0xfc,
                    0xf5, 0xe0, 0x95, 0xd5, 0xd0, 0x0a, 0xf6, 0x00, 0xdb, 0x18, 0xcb, 0x2c, 0x04,
                    0xb3, 0xed, 0xd0, 0x3c, 0xc7, 0x44, 0xa2, 0x88, 0x8a, 0xe4, 0x0c, 0xaa, 0x23,
                    0x29, 0x46, 0xc5, 0xe7, 0xe1,
                ],
            ),
            (
                12345,
                [
                    0x05, 0x30, 0xc1, 0xbd, 0xc4, 0xcd, 0x6b, 0x14, 0x08, 0xbe, 0x09, 0x33, 0xc4,
                    0xa4, 0x1a, 0xc3, 0x51, 0x33, 0x50, 0xee, 0xf3, 0x68, 0x50, 0xb8, 0x04, 0x70,
                    0x8e, 0x1f, 0x33, 0x89, 0x32, 0xce, 0x01, 0xb6, 0x55, 0xa1, 0x63, 0x34, 0x4a,
                    0x45, 0x00, 0xb2, 0x81, 0xc8, 0x75, 0x0c, 0x46, 0x1f, 0x00, 0x38, 0xe7, 0x6f,
                    0x31, 0xb5, 0xae, 0xf9, 0xd7, 0xc8, 0xf1, 0x61, 0x6d, 0x24, 0x46, 0xfb, 0x1f,
                    0x03, 0x80, 0xac, 0xa5, 0x86, 0xb9, 0x26, 0x8a, 0x11, 0x5e, 0x8d, 0x19, 0x1c,
                    0x5e, 0x47, 0xed, 0x04, 0xc4, 0xb7, 0x7b, 0x72, 0x39, 0x47, 0x40, 0x02, 0x57,
                    0x06, 0x84, 0x5c, 0x9c, 0xb7,
                ],
            ),
        ];

        for (k, expected) in COMPRESSED {
            let p = (Point::GENERATOR * &s(*k)).to_affine().unwrap();
            assert_eq!(&p.to_compressed(), expected, "k = {}", k);
            assert_eq!(PointAffine::from_compressed(expected).unwrap(), p);
        }
        for (k, expected) in UNCOMPRESSED {
            let p = (Point::GENERATOR * &s(*k)).to_affine().unwrap();
            assert_eq!(&p.to_uncompressed(), expected, "k = {}", k);
            assert_eq!(PointAffine::from_uncompressed(expected).unwrap(), p);
        }
    }
}
