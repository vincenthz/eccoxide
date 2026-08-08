//! BLS12-381 G2
//!
//! G2 is the prime-order-`r` subgroup of the sextic twist
//! `E'(Fp2): y^2 = x^3 + 4(1 + u)` (so `a = 0`, `b = 4 + 4u`) defined over the
//! quadratic extension field [`Fp2`](super::fp2::Fp2). Its scalar field is the
//! same [`Scalar`](super::scalar::Scalar) as G1.
//!
//! Like [`super::g1`], this wires the generic curve machinery (affine /
//! projective points, complete addition and the fixed-base comb) into the
//! [`CurveGroup`] abstraction, this time over `Fp2`.
//!
//! # G2 versus the twist
//!
//! `E'(Fp2)` has order `h₂·r` with a 507-bit cofactor `h₂`, so — even more so
//! than in [`super::g1`] — satisfying the twist equation is much weaker than
//! being in the prime-order group G2 the pairing is defined on.
//! [`Point::is_in_subgroup`] tests the difference, and
//! [`Point::clear_cofactor`] establishes it.
//!
//! [`CurveGroup`]: crate::curve::group::CurveGroup

use crate::curve::bls12_381::fp2::Fp2 as FieldElement;
#[cfg(feature = "bls12-381-hash-to-curve")]
use crate::curve::bls12_381::hash_to_curve;
use crate::curve::bls12_381::scalar::Scalar;
use crate::curve::bls12_381::BLS_X;
use crate::curve::field::Sign;
use crate::curve::{
    affine, projective,
    weierstrass::{WeierstrassCurve, WeierstrassCurveA0},
};
use crate::mp::ct::{Choice, CtEqual};
use crate::params::bls12_381::g2::*;
use crate::params::bls12_381::ORDER_BYTES;
use crate::{
    bls12_381_define_point_serialization, fiat_define_weierstrass_curve,
    fiat_define_weierstrass_curve_a0, fiat_define_weierstrass_points,
};

fiat_define_weierstrass_curve!(FieldElement);
fiat_define_weierstrass_curve_a0!(FieldElement);
fiat_define_weierstrass_points!(FieldElement);
bls12_381_define_point_serialization!(FieldElement);

/// `ξ^-((p-1)/3)`, the x-coordinate coefficient of `ψ`; see [`PSI_X_COEFF_BYTES`].
const PSI_X_COEFF: FieldElement = FieldElement::from_bytes_unchecked(&PSI_X_COEFF_BYTES);

/// `ξ^-((p-1)/2)`, the y-coordinate coefficient of `ψ`; see [`PSI_Y_COEFF_BYTES`].
const PSI_Y_COEFF: FieldElement = FieldElement::from_bytes_unchecked(&PSI_Y_COEFF_BYTES);

impl Point {
    /// The untwist-Frobenius-twist endomorphism `ψ` of the twist.
    ///
    /// `ψ = φ ∘ π ∘ φ⁻¹`, where `π` is the `p`-power Frobenius of `E(Fp12)` and
    /// `φ(x, y) = (x·w², y·w³)` is the twist map of
    /// [`super::pairing`](super::pairing). Both `w` powers collapse to `Fp2`
    /// constants, `w^(2(p-1)) = ξ^((p-1)/3)` and `w^(3(p-1)) = ξ^((p-1)/2)` for
    /// the tower non-residue `ξ = 1 + u`, which leaves
    ///
    /// ```text
    /// ψ(x, y) = (x^p · ξ^-((p-1)/3), y^p · ξ^-((p-1)/2))
    /// ```
    ///
    /// i.e. a conjugation and one `Fp2` multiplication per coordinate — no
    /// `Fp12` arithmetic despite the detour through it. On G2, `ψ` is
    /// multiplication by `p`, and `p ≡ x (mod r)`, which is what
    /// [`Self::is_in_subgroup`] exploits.
    fn psi(&self) -> Self {
        Point(projective::Point {
            x: &self.0.x.frobenius_map() * &PSI_X_COEFF,
            y: &self.0.y.frobenius_map() * &PSI_Y_COEFF,
            z: self.0.z.frobenius_map(),
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

    /// Whether the point lies in the prime-order-`r` subgroup G2, and not merely
    /// on the twist.
    ///
    /// Points read from an untrusted source have to be checked: the cofactor of
    /// `E'(Fp2)` is 507 bits wide, so all but a vanishing fraction of the twist's
    /// points are outside G2, and feeding one to the pairing breaks the
    /// assumptions the protocols above it rest on. The standard decoders
    /// ([`PointAffine::from_compressed`] and friends) already run this check.
    ///
    /// The test is the one of section 4 of <https://eprint.iacr.org/2021/1130>
    /// (proof of correctness revised in <https://eprint.iacr.org/2022/352>): for
    /// `Q` on the twist,
    ///
    /// ```text
    /// Q ∈ G2   ⟺   ψ(Q) = [x]Q
    /// ```
    ///
    /// `ψ` being the untwist-Frobenius-twist endomorphism
    /// `ψ(x, y) = (x^p·ξ^-((p-1)/3), y^p·ξ^-((p-1)/2))`, which acts on G2 as
    /// `[p]`, and `p ≡ x (mod r)`. Since the seed is only 64 bits, this costs one
    /// multiplication by it plus the two `Fp2` multiplications of `ψ`, about 5x
    /// less than the multiplication by the 255-bit `r` the definition
    /// `[r]Q == 𝒪` would ask for.
    ///
    /// The identity is in the subgroup, so this is true of it.
    pub fn is_in_subgroup(&self) -> Choice {
        // [x]Q = -[|x|]Q, the seed being negative
        let xq = -&self.mul_by_abs_x();
        self.psi().0.ct_eq(&xq.0)
    }

    /// `ψ` applied twice.
    ///
    /// The two Frobenius applications collapse (`a^(p²) = a` in `Fp2`), so this
    /// is really multiplication of each coordinate by the norm of the
    /// corresponding `ψ` coefficient; it is spelled as two `ψ` here since the
    /// cofactor clearing it serves is dominated by the seed multiplications.
    fn psi2(&self) -> Self {
        self.psi().psi()
    }

    /// Send an arbitrary point of the twist into the prime-order subgroup G2.
    ///
    /// This is the last step of hash-to-curve: a map onto `E'(Fp2)` lands
    /// anywhere on the twist, and this brings the result into the group the
    /// pairing is defined on. The output is in G2 for *any* input, the identity
    /// included, and the map is surjective onto G2.
    ///
    /// The cofactor of the twist is 507 bits wide, so multiplying by it directly
    /// would cost far more than the group operations around it. Instead this is
    /// the endomorphism chain of <https://eprint.iacr.org/2017/419>,
    ///
    /// ```text
    /// [x² - x - 1]Q + [x - 1]ψ(Q) + [2]ψ²(Q)
    /// ```
    ///
    /// evaluated as `ψ²([2]Q) + [x]([x]Q + ψ(Q)) - [x]Q - ψ(Q) - Q`, i.e. two
    /// multiplications by the 64-bit seed and three `ψ` applications. It works
    /// out to multiplication by `3(x² - 1)·h₂`, a multiple of the cofactor `h₂`
    /// rather than `h₂` itself — which clears it just as well, and is exactly
    /// the `h_eff` of section 8.8.2 of RFC 9380.
    ///
    /// Note this is not the identity on points that are already in G2: there `ψ`
    /// is `[x]`, so it is multiplication by the scalar `4x² - 2x - 1`. Use
    /// [`Self::is_in_subgroup`] to *test* membership; this one *establishes* it.
    ///
    /// The result is projective because clearing the cofactor of a point of
    /// order dividing `h₂` gives the identity, which the affine type cannot hold.
    pub fn clear_cofactor(&self) -> Self {
        let xq = -&self.mul_by_abs_x(); // [x]Q
        let psi_q = self.psi(); // ψ(Q)
        let x_sum = -&(&xq + &psi_q).mul_by_abs_x(); // [x²]Q + [x]ψ(Q)
        let psi2_2q = Point(self.0.double_a0::<Curve>()).psi2(); // ψ²([2]Q)

        let acc = &psi2_2q + &x_sum;
        let acc = &acc - &xq;
        let acc = &acc - &psi_q;
        &acc - self
    }
}

impl PointAffine {
    /// Whether the point lies in the prime-order-`r` subgroup G2, and not merely
    /// on the twist; see [`Point::is_in_subgroup`], which this defers to.
    pub fn is_in_subgroup(&self) -> Choice {
        Point::from(self).is_in_subgroup()
    }

    /// Send an arbitrary point of the twist into the prime-order subgroup G2;
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

    /// Hash a message to a point of G2, as the RFC 9380 suite
    /// `BLS12381G2_XMD:SHA-256_SSWU_RO_`.
    ///
    /// The result is uniformly distributed over G2 and the map is
    /// indifferentiable from a random oracle, which is what protocols building
    /// on it (BLS signatures, verifiable random functions, commitments to
    /// arbitrary data) assume. Prefer this over
    /// [`Self::encode_to_curve`] unless the input is already uniform.
    ///
    /// `dst` is the domain separation tag: pick one string per application and
    /// per use within it, and keep it fixed. RFC 9380 section 3.1 recommends
    /// embedding the suite name, e.g.
    /// `b"MYPROTOCOL-V01-CS01-with-BLS12381G2_XMD:SHA-256_SSWU_RO_"`. A tag longer
    /// than 255 bytes is hashed down as the specification prescribes, so any
    /// length works.
    ///
    /// The message is hashed to two field elements, each mapped to the twist
    /// and added, and the sum has its cofactor cleared — see the
    /// [`hash_to_curve` module](super::hash_to_curve) for the pipeline.
    pub fn hash_to_curve(msg: &[u8], dst: &[u8]) -> Self {
        let [u0, u1] = hash_to_curve::hash_to_field_g2(msg, dst);
        let q0 = Point(hash_to_curve::map_to_curve_g2(&u0));
        let q1 = Point(hash_to_curve::map_to_curve_g2(&u1));
        (&q0 + &q1).clear_cofactor()
    }

    /// Encode a message to a point of G2, as the RFC 9380 suite
    /// `BLS12381G2_XMD:SHA-256_SSWU_NU_`.
    ///
    /// Cheaper than [`Self::hash_to_curve`] — one field element and one map
    /// instead of two — but *non-uniform*: its image is a fraction of G2,
    /// and it is not indifferentiable from a random oracle. Only use it when
    /// the caller already guarantees a uniform input, or when the protocol
    /// explicitly asks for the `NU_` suite.
    ///
    /// See [`Self::hash_to_curve`] about the domain separation tag.
    pub fn encode_to_curve(msg: &[u8], dst: &[u8]) -> Self {
        let u = hash_to_curve::encode_to_field_g2(msg, dst);
        Point(hash_to_curve::map_to_curve_g2(&u)).clear_cofactor()
    }
}

#[cfg(test)]
mod tests {
    use super::{Point, PointAffine, Scalar};
    use crate::curve::group::CurveGroup;

    fn s(v: u64) -> Scalar {
        Scalar::from_u64(v)
    }

    #[test]
    fn generator_on_curve() {
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
        // the G2 fixed-base comb table must agree with variable-base multiplication
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

    crate::bls12_381_point_serialization_unittest!();

    /// Prime-order-subgroup membership, and the endomorphism it is built on.
    mod subgroup {
        use super::super::{Curve, FieldElement, Point, PointAffine};
        use super::{s, Scalar};
        use crate::curve::bls12_381::BLS_X;
        use crate::curve::field::Sign;
        use crate::params::bls12_381::ORDER_BYTES;

        /// Points of the twist that are **not** in the prime-order subgroup, as their
        /// `(compressed, uncompressed)` encodings.
        ///
        /// The first two are the decompressions of the first two x-coordinates
        /// of the form `i + (i + 1)u` that are on the twist at all, `x = 1 + 2u`
        /// and `x = 4 + 5u`, with opposite y-signs so that the sort flag is
        /// pinned both ways; the third is `G` plus the first, the shape an
        /// attacker sends — a genuine G2 point with a torsion component added.
        ///
        /// Generated, and checked against `[r]Q != 𝒪`, by `sage/bls12_381.sage`.
        const OFF_SUBGROUP: &[([u8; 96], [u8; 192])] = &[
            (
                [
                    0xa0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x01,
                ],
                [
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x01, 0x0f, 0x9a, 0x36, 0xa7, 0x74, 0xaa, 0x58, 0x43,
                    0x2e, 0x4d, 0x62, 0x04, 0x53, 0xde, 0xa1, 0xa4, 0xbc, 0x03, 0xb7, 0x99, 0xd0,
                    0xdb, 0x1d, 0x4f, 0x82, 0x0b, 0xbd, 0xfc, 0x17, 0xa6, 0x37, 0x6f, 0x80, 0xed,
                    0xb7, 0x77, 0xdb, 0x21, 0x52, 0xa0, 0xeb, 0x5c, 0x3e, 0x87, 0xf0, 0x23, 0x55,
                    0xac, 0x01, 0x67, 0x46, 0xf5, 0xf8, 0x1e, 0x24, 0x89, 0xa9, 0x56, 0xb7, 0x2c,
                    0xeb, 0xb4, 0xa2, 0x3e, 0x9a, 0xd2, 0xf3, 0x37, 0x48, 0x8b, 0xc9, 0x54, 0x0c,
                    0x19, 0x0c, 0x6b, 0x89, 0x45, 0x14, 0x34, 0x85, 0x8a, 0x81, 0x97, 0x7d, 0xff,
                    0xa1, 0x28, 0x6d, 0x83, 0xc5, 0x7d, 0xf5, 0xf9, 0x2a, 0xe7,
                ],
            ),
            (
                [
                    0x80, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x05, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x04,
                ],
                [
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x05, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x04, 0x02, 0x04, 0xd9, 0xad, 0xdf, 0x67, 0xf0, 0xfa,
                    0xa8, 0x94, 0xfb, 0xf8, 0x39, 0xbf, 0xd2, 0x69, 0xe0, 0x52, 0x10, 0x88, 0x74,
                    0xcc, 0x4b, 0x13, 0x35, 0x58, 0xed, 0xea, 0xcc, 0xb0, 0x80, 0x16, 0x63, 0x83,
                    0xfa, 0xf7, 0xa2, 0x70, 0x85, 0x09, 0x05, 0x33, 0x0b, 0xd9, 0x58, 0xee, 0x4e,
                    0x67, 0x0d, 0x2e, 0x09, 0x08, 0x3c, 0x74, 0x76, 0xd3, 0x42, 0xc8, 0x30, 0xe2,
                    0xc4, 0x3a, 0x6e, 0x52, 0x2a, 0x87, 0xf3, 0x08, 0x3c, 0xd5, 0xfd, 0x6c, 0x34,
                    0x2b, 0xfd, 0x46, 0x99, 0x48, 0xf4, 0x39, 0xe8, 0x83, 0xa0, 0xe3, 0x26, 0xfd,
                    0xa6, 0xa6, 0x79, 0x9e, 0x4d, 0xfc, 0x8e, 0x64, 0xde, 0x87,
                ],
            ),
            (
                [
                    0xab, 0x00, 0xe0, 0xea, 0x1f, 0x49, 0xcb, 0x65, 0x5d, 0x95, 0x4d, 0xd5, 0x18,
                    0x53, 0xc0, 0x39, 0x17, 0xb9, 0xc1, 0xe3, 0xd4, 0xcb, 0x68, 0xcf, 0x55, 0x65,
                    0x48, 0x26, 0x3e, 0x34, 0xd3, 0x3b, 0xb9, 0x05, 0x3b, 0x44, 0xfe, 0x7c, 0x12,
                    0x6f, 0x35, 0x4e, 0x83, 0xd2, 0xb9, 0xd7, 0x16, 0xae, 0x08, 0x74, 0x25, 0x63,
                    0xd7, 0x98, 0x66, 0xc4, 0x3b, 0xe0, 0x4f, 0xec, 0x10, 0x65, 0xe2, 0x49, 0xde,
                    0xcc, 0xf5, 0xd0, 0xcf, 0x93, 0x14, 0xa8, 0xb5, 0xf1, 0xe3, 0x20, 0xb4, 0x10,
                    0xe1, 0xb9, 0xd7, 0x52, 0xbf, 0xdd, 0xa7, 0xd5, 0x9a, 0x89, 0x74, 0xe6, 0xfd,
                    0x99, 0x85, 0x2a, 0x4b, 0x20,
                ],
                [
                    0x0b, 0x00, 0xe0, 0xea, 0x1f, 0x49, 0xcb, 0x65, 0x5d, 0x95, 0x4d, 0xd5, 0x18,
                    0x53, 0xc0, 0x39, 0x17, 0xb9, 0xc1, 0xe3, 0xd4, 0xcb, 0x68, 0xcf, 0x55, 0x65,
                    0x48, 0x26, 0x3e, 0x34, 0xd3, 0x3b, 0xb9, 0x05, 0x3b, 0x44, 0xfe, 0x7c, 0x12,
                    0x6f, 0x35, 0x4e, 0x83, 0xd2, 0xb9, 0xd7, 0x16, 0xae, 0x08, 0x74, 0x25, 0x63,
                    0xd7, 0x98, 0x66, 0xc4, 0x3b, 0xe0, 0x4f, 0xec, 0x10, 0x65, 0xe2, 0x49, 0xde,
                    0xcc, 0xf5, 0xd0, 0xcf, 0x93, 0x14, 0xa8, 0xb5, 0xf1, 0xe3, 0x20, 0xb4, 0x10,
                    0xe1, 0xb9, 0xd7, 0x52, 0xbf, 0xdd, 0xa7, 0xd5, 0x9a, 0x89, 0x74, 0xe6, 0xfd,
                    0x99, 0x85, 0x2a, 0x4b, 0x20, 0x16, 0x30, 0xad, 0x86, 0x89, 0x4d, 0x6e, 0x8c,
                    0xe9, 0x0b, 0xdd, 0xd6, 0xd2, 0xc5, 0xc2, 0x9f, 0x3a, 0xc3, 0xe2, 0x2c, 0x2b,
                    0x67, 0x6b, 0xde, 0xca, 0x72, 0x55, 0x0f, 0xcd, 0x80, 0x77, 0xca, 0x9c, 0x68,
                    0x11, 0x1c, 0x12, 0xce, 0xe6, 0x8e, 0x1f, 0xa2, 0xbf, 0x67, 0x88, 0xc6, 0xf8,
                    0x3e, 0x15, 0xdb, 0x9c, 0xf6, 0x99, 0xaf, 0xb1, 0xe3, 0x23, 0x3b, 0x38, 0x66,
                    0x80, 0xf1, 0x35, 0x45, 0x63, 0x11, 0xd3, 0xec, 0xf4, 0xf4, 0x84, 0x6b, 0xac,
                    0xc2, 0x56, 0x28, 0xb7, 0xe2, 0x6c, 0xc5, 0x44, 0x91, 0x66, 0x5e, 0x9c, 0x32,
                    0xd5, 0x40, 0x06, 0xdb, 0xfd, 0x72, 0x78, 0x1d, 0xb6, 0x5e,
                ],
            ),
        ];

        /// `[r]q`, the definition the fast check has to agree with: a point of
        /// the twist is in the prime-order subgroup exactly when the group order
        /// kills it. Deliberately computed the slow way, by a generic scalar
        /// multiplication by `r` — which no `Scalar` could hold, `r` being the
        /// modulus of that very field.
        fn mul_by_order(q: &Point) -> Point {
            Point(q.0.scale_a0::<Curve>(&ORDER_BYTES))
        }

        /// Points of the twist found by decompressing small x-coordinates (real
        /// ones, `c1 = 0`). The cofactor of `E'(Fp2)` is 380 bits wide, so none
        /// of these is in G2 by any plausible chance, and about half the x tried
        /// are on the twist.
        fn small_x_points() -> Vec<PointAffine> {
            let points: Vec<_> = (1..32u64)
                .filter_map(|i| {
                    PointAffine::decompress(&FieldElement::from(i), Sign::Positive).into_option()
                })
                .collect();
            assert!(!points.is_empty(), "no small x gave a point of the twist");
            points
        }

        /// Multiples of the generator are in G2, whichever representation they
        /// are held in, and so is the identity.
        #[test]
        fn generator_multiples_are_in_the_subgroup() {
            for v in [1u64, 2, 3, 12345, 0xdead_beef] {
                let q = Point::GENERATOR * &s(v);
                assert!(q.is_in_subgroup().is_true(), "[{}]G rejected", v);
                assert!(q.to_affine().unwrap().is_in_subgroup().is_true());
            }
            assert!(Point::INFINITY.is_in_subgroup().is_true());
        }

        /// The check must agree with `[r]Q == 𝒪` on points of the twist, both in
        /// the subgroup and outside it — that equivalence is the whole claim of
        /// the `ψ` test.
        #[test]
        fn matches_multiplication_by_the_order() {
            let mut points: Vec<Point> = [1u64, 3, 12345]
                .iter()
                .map(|v| Point::GENERATOR * &s(*v))
                .collect();
            points.push(Point::INFINITY);
            points.extend(small_x_points().iter().map(Point::from));
            // the shape an attacker would send: a point of G2 with a torsion
            // component added, so the order is a proper multiple of r
            points.extend(
                small_x_points()
                    .iter()
                    .map(|q| &Point::GENERATOR + &Point::from(q)),
            );

            let (mut inside, mut outside) = (0, 0);
            for q in points {
                let expected = mul_by_order(&q) == Point::INFINITY;
                assert_eq!(
                    q.is_in_subgroup().is_true(),
                    expected,
                    "the psi check disagrees with [r]Q == O"
                );
                if expected {
                    inside += 1;
                } else {
                    outside += 1;
                }
            }
            assert!(inside > 0 && outside > 0, "one side was never exercised");
        }

        /// The eigenvalue of `ψ` on G2, which is what pins its two coefficients:
        /// `ψ` acts as `[p]`, and `p ≡ x (mod r)` for the (negative) seed.
        #[test]
        fn psi_is_multiplication_by_the_seed() {
            let x = -&Scalar::from_u64(BLS_X);
            for v in [1u64, 2, 12345] {
                let q = Point::GENERATOR * &s(v);
                assert_eq!(q.psi(), &q * &x, "psi != [x] on [{}]G", v);
            }
        }

        /// `ψ` is an endomorphism of the whole twist, not just of G2.
        #[test]
        fn psi_stays_on_the_twist() {
            for q in small_x_points() {
                let r = Point::from(&q).psi().to_affine().unwrap();
                let (x, y) = r.to_coordinate();
                assert!(PointAffine::from_coordinate(x, y).is_some());
            }
        }

        /// The double-and-add over the seed agrees with a generic scalar
        /// multiplication by it.
        #[test]
        fn mul_by_abs_x_matches_scalar_mul() {
            let k = Scalar::from_u64(BLS_X);
            for v in [1u64, 7, 12345] {
                let q = Point::GENERATOR * &s(v);
                assert_eq!(q.mul_by_abs_x(), &q * &k);
            }
        }

        /// `h_eff` for G2, `3(x² - 1)·h₂`, as section 8.8.2 of RFC 9380 gives it
        /// (big-endian). The endomorphism chain of [`Point::clear_cofactor`] has
        /// to agree with a plain multiplication by this.
        const H_EFF: [u8; 80] = [
            0x0b, 0xc6, 0x9f, 0x08, 0xf2, 0xee, 0x75, 0xb3, 0x58, 0x4c, 0x6a, 0x0e, 0xa9, 0x1b,
            0x35, 0x28, 0x88, 0xe2, 0xa8, 0xe9, 0x14, 0x5a, 0xd7, 0x68, 0x99, 0x86, 0xff, 0x03,
            0x15, 0x08, 0xff, 0xe1, 0x32, 0x9c, 0x2f, 0x17, 0x87, 0x31, 0xdb, 0x95, 0x6d, 0x82,
            0xbf, 0x01, 0x5d, 0x12, 0x12, 0xb0, 0x2e, 0xc0, 0xec, 0x69, 0xd7, 0x47, 0x7c, 0x1a,
            0xe9, 0x54, 0xcb, 0xc0, 0x66, 0x89, 0xf6, 0xa3, 0x59, 0x89, 0x4c, 0x0a, 0xde, 0xbb,
            0xf6, 0xb4, 0xe8, 0x02, 0x00, 0x05, 0xaa, 0xa9, 0x55, 0x51,
        ];

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

            for q in points {
                let cleared = q.clear_cofactor();
                assert!(cleared.is_in_subgroup().is_true(), "cleared point rejected");
                assert_eq!(
                    mul_by_order(&cleared),
                    Point::INFINITY,
                    "[r] does not kill the cleared point"
                );
            }
        }

        /// It is multiplication by the `h_eff` of RFC 9380 section 8.8.2,
        /// `3(x² - 1)·h₂`, and not by anything else — checked against a generic
        /// multiplication by those bytes, for points on and off the subgroup.
        #[test]
        fn clear_cofactor_is_multiplication_by_h_eff() {
            let mut points: Vec<Point> = vec![Point::GENERATOR, Point::GENERATOR * &s(12345)];
            points.extend(small_x_points().iter().map(Point::from));

            for q in points {
                let expected = Point(q.0.scale_a0::<Curve>(&H_EFF));
                assert_eq!(q.clear_cofactor(), expected);
            }
            // on the subgroup psi is [x], so the chain is the scalar 4x^2-2x-1
            // there: not a no-op, and congruent to h_eff modulo r
            let x = -&Scalar::from_u64(BLS_X);
            let four = Scalar::from_u64(4);
            let two = Scalar::from_u64(2);
            let lambda = &(&(&four * &(&x * &x)) - &(&two * &x)) - &Scalar::from_u64(1);
            let g = Point::GENERATOR;
            assert_eq!(g.clear_cofactor(), &g * &lambda);
            assert_ne!(g.clear_cofactor(), g);
        }

        /// The identity clears to the identity, a stray point does not collapse
        /// to it, and the affine entry point agrees with the projective one.
        #[test]
        fn clear_cofactor_edge_cases() {
            assert_eq!(Point::INFINITY.clear_cofactor(), Point::INFINITY);
            for q in small_x_points() {
                assert_ne!(
                    Point::from(&q).clear_cofactor(),
                    Point::INFINITY,
                    "a twist point cleared to the identity"
                );
                assert_eq!(q.clear_cofactor(), Point::from(&q).clear_cofactor());
            }
        }
    }

    /// Known-answer vectors for the standard G2 encodings: `k * G` for small
    /// `k`, byte for byte.
    ///
    /// Note the `c1 || c0` order of each Fp2 coordinate: the `k = 1` vector,
    /// the published encoding of the generator, starts with its
    /// `x.c1 = 0x13e02b60...` and not with its `x.c0 = 0x024aa2b2...`. Its sort
    /// flag is clear, while `k = 2` and `k = 0xdeadbeef` have it set, so the
    /// vectors pin the flag both ways. The multiples were computed
    /// independently.
    #[test]
    fn serialization_kat() {
        const COMPRESSED: &[(u64, [u8; 96])] = &[
            (
                1,
                [
                    0x93, 0xe0, 0x2b, 0x60, 0x52, 0x71, 0x9f, 0x60, 0x7d, 0xac, 0xd3, 0xa0, 0x88,
                    0x27, 0x4f, 0x65, 0x59, 0x6b, 0xd0, 0xd0, 0x99, 0x20, 0xb6, 0x1a, 0xb5, 0xda,
                    0x61, 0xbb, 0xdc, 0x7f, 0x50, 0x49, 0x33, 0x4c, 0xf1, 0x12, 0x13, 0x94, 0x5d,
                    0x57, 0xe5, 0xac, 0x7d, 0x05, 0x5d, 0x04, 0x2b, 0x7e, 0x02, 0x4a, 0xa2, 0xb2,
                    0xf0, 0x8f, 0x0a, 0x91, 0x26, 0x08, 0x05, 0x27, 0x2d, 0xc5, 0x10, 0x51, 0xc6,
                    0xe4, 0x7a, 0xd4, 0xfa, 0x40, 0x3b, 0x02, 0xb4, 0x51, 0x0b, 0x64, 0x7a, 0xe3,
                    0xd1, 0x77, 0x0b, 0xac, 0x03, 0x26, 0xa8, 0x05, 0xbb, 0xef, 0xd4, 0x80, 0x56,
                    0xc8, 0xc1, 0x21, 0xbd, 0xb8,
                ],
            ),
            (
                2,
                [
                    0xaa, 0x4e, 0xde, 0xf9, 0xc1, 0xed, 0x7f, 0x72, 0x9f, 0x52, 0x0e, 0x47, 0x73,
                    0x0a, 0x12, 0x4f, 0xd7, 0x06, 0x62, 0xa9, 0x04, 0xba, 0x10, 0x74, 0x72, 0x81,
                    0x14, 0xd1, 0x03, 0x1e, 0x15, 0x72, 0xc6, 0xc8, 0x86, 0xf6, 0xb5, 0x7e, 0xc7,
                    0x2a, 0x61, 0x78, 0x28, 0x8c, 0x47, 0xc3, 0x35, 0x77, 0x16, 0x38, 0x53, 0x39,
                    0x57, 0xd5, 0x40, 0xa9, 0xd2, 0x37, 0x0f, 0x17, 0xcc, 0x7e, 0xd5, 0x86, 0x3b,
                    0xc0, 0xb9, 0x95, 0xb8, 0x82, 0x5e, 0x0e, 0xe1, 0xea, 0x1e, 0x1e, 0x4d, 0x00,
                    0xdb, 0xae, 0x81, 0xf1, 0x4b, 0x0b, 0xf3, 0x61, 0x1b, 0x78, 0xc9, 0x52, 0xaa,
                    0xca, 0xb8, 0x27, 0xa0, 0x53,
                ],
            ),
            (
                3,
                [
                    0x89, 0x38, 0x02, 0x75, 0xbb, 0xc8, 0xe5, 0xdc, 0xea, 0x7d, 0xc4, 0xdd, 0x7e,
                    0x05, 0x50, 0xff, 0x2a, 0xc4, 0x80, 0x90, 0x53, 0x96, 0xed, 0xa5, 0x50, 0x62,
                    0x65, 0x0f, 0x8d, 0x25, 0x1c, 0x96, 0xeb, 0x48, 0x06, 0x73, 0x93, 0x7c, 0xc6,
                    0xd9, 0xd6, 0xa4, 0x4a, 0xaa, 0x56, 0xca, 0x66, 0xdc, 0x12, 0x29, 0x15, 0xc8,
                    0x24, 0xa0, 0x85, 0x7e, 0x2e, 0xe4, 0x14, 0xa3, 0xdc, 0xcb, 0x23, 0xae, 0x69,
                    0x1a, 0xe5, 0x43, 0x29, 0x78, 0x13, 0x15, 0xa0, 0xc7, 0x5d, 0xf1, 0xc0, 0x4d,
                    0x6d, 0x7a, 0x50, 0xa0, 0x30, 0xfc, 0x86, 0x6f, 0x09, 0xd5, 0x16, 0x02, 0x0e,
                    0xf8, 0x23, 0x24, 0xaf, 0xae,
                ],
            ),
            (
                12345,
                [
                    0x84, 0x9d, 0x5b, 0x3d, 0x40, 0xfe, 0x47, 0x5b, 0x14, 0x5e, 0xeb, 0xf5, 0x3d,
                    0x97, 0x98, 0x1b, 0xde, 0x5a, 0x64, 0xde, 0xa2, 0x96, 0x48, 0x07, 0xf8, 0x25,
                    0x61, 0xe7, 0x09, 0xe8, 0x04, 0xfe, 0xe3, 0xec, 0xfb, 0x53, 0x56, 0x63, 0x1b,
                    0x2d, 0xed, 0xbe, 0x82, 0xd3, 0xd1, 0xda, 0xd0, 0xbb, 0x03, 0x7e, 0xce, 0x3e,
                    0xcc, 0x51, 0x22, 0x26, 0xa1, 0xe5, 0x6f, 0xbe, 0x0b, 0x33, 0xaa, 0xb2, 0x08,
                    0x0a, 0xb4, 0x67, 0xd1, 0x4a, 0xad, 0xef, 0xf5, 0xdc, 0xd8, 0xad, 0xc6, 0x61,
                    0x3b, 0x92, 0x6b, 0xc9, 0x76, 0x01, 0xa4, 0xa1, 0xf1, 0x28, 0x77, 0x93, 0x75,
                    0x7b, 0x10, 0xd6, 0x8a, 0x93,
                ],
            ),
            (
                0xdead_beef,
                [
                    0xa2, 0x91, 0x08, 0x66, 0xed, 0xe4, 0x27, 0x89, 0xc6, 0x85, 0xaf, 0x18, 0x62,
                    0x34, 0xbf, 0xe4, 0x66, 0x4e, 0x27, 0x0e, 0xd4, 0x47, 0xd2, 0xd1, 0x42, 0xa1,
                    0x59, 0xcf, 0x0c, 0xf9, 0xba, 0xc6, 0xc3, 0x78, 0x9d, 0x55, 0x8e, 0xa2, 0xe2,
                    0x4b, 0x60, 0xe4, 0x8c, 0xeb, 0x8b, 0x96, 0x0f, 0xa8, 0x11, 0x5f, 0xd9, 0x14,
                    0x42, 0x73, 0x34, 0xb6, 0xc4, 0x5a, 0x82, 0x3e, 0x1e, 0x01, 0x1b, 0x6b, 0x2a,
                    0xa5, 0x13, 0xf9, 0x5f, 0xd6, 0x4c, 0x95, 0x61, 0x28, 0xfd, 0x63, 0x18, 0x57,
                    0x94, 0x5a, 0xdf, 0x8a, 0xa7, 0xb0, 0x36, 0x4c, 0xee, 0xa7, 0x37, 0xc0, 0x67,
                    0x1f, 0xd9, 0x39, 0x21, 0x9e,
                ],
            ),
        ];
        const UNCOMPRESSED: &[(u64, [u8; 192])] = &[
            (
                1,
                [
                    0x13, 0xe0, 0x2b, 0x60, 0x52, 0x71, 0x9f, 0x60, 0x7d, 0xac, 0xd3, 0xa0, 0x88,
                    0x27, 0x4f, 0x65, 0x59, 0x6b, 0xd0, 0xd0, 0x99, 0x20, 0xb6, 0x1a, 0xb5, 0xda,
                    0x61, 0xbb, 0xdc, 0x7f, 0x50, 0x49, 0x33, 0x4c, 0xf1, 0x12, 0x13, 0x94, 0x5d,
                    0x57, 0xe5, 0xac, 0x7d, 0x05, 0x5d, 0x04, 0x2b, 0x7e, 0x02, 0x4a, 0xa2, 0xb2,
                    0xf0, 0x8f, 0x0a, 0x91, 0x26, 0x08, 0x05, 0x27, 0x2d, 0xc5, 0x10, 0x51, 0xc6,
                    0xe4, 0x7a, 0xd4, 0xfa, 0x40, 0x3b, 0x02, 0xb4, 0x51, 0x0b, 0x64, 0x7a, 0xe3,
                    0xd1, 0x77, 0x0b, 0xac, 0x03, 0x26, 0xa8, 0x05, 0xbb, 0xef, 0xd4, 0x80, 0x56,
                    0xc8, 0xc1, 0x21, 0xbd, 0xb8, 0x06, 0x06, 0xc4, 0xa0, 0x2e, 0xa7, 0x34, 0xcc,
                    0x32, 0xac, 0xd2, 0xb0, 0x2b, 0xc2, 0x8b, 0x99, 0xcb, 0x3e, 0x28, 0x7e, 0x85,
                    0xa7, 0x63, 0xaf, 0x26, 0x74, 0x92, 0xab, 0x57, 0x2e, 0x99, 0xab, 0x3f, 0x37,
                    0x0d, 0x27, 0x5c, 0xec, 0x1d, 0xa1, 0xaa, 0xa9, 0x07, 0x5f, 0xf0, 0x5f, 0x79,
                    0xbe, 0x0c, 0xe5, 0xd5, 0x27, 0x72, 0x7d, 0x6e, 0x11, 0x8c, 0xc9, 0xcd, 0xc6,
                    0xda, 0x2e, 0x35, 0x1a, 0xad, 0xfd, 0x9b, 0xaa, 0x8c, 0xbd, 0xd3, 0xa7, 0x6d,
                    0x42, 0x9a, 0x69, 0x51, 0x60, 0xd1, 0x2c, 0x92, 0x3a, 0xc9, 0xcc, 0x3b, 0xac,
                    0xa2, 0x89, 0xe1, 0x93, 0x54, 0x86, 0x08, 0xb8, 0x28, 0x01,
                ],
            ),
            (
                12345,
                [
                    0x04, 0x9d, 0x5b, 0x3d, 0x40, 0xfe, 0x47, 0x5b, 0x14, 0x5e, 0xeb, 0xf5, 0x3d,
                    0x97, 0x98, 0x1b, 0xde, 0x5a, 0x64, 0xde, 0xa2, 0x96, 0x48, 0x07, 0xf8, 0x25,
                    0x61, 0xe7, 0x09, 0xe8, 0x04, 0xfe, 0xe3, 0xec, 0xfb, 0x53, 0x56, 0x63, 0x1b,
                    0x2d, 0xed, 0xbe, 0x82, 0xd3, 0xd1, 0xda, 0xd0, 0xbb, 0x03, 0x7e, 0xce, 0x3e,
                    0xcc, 0x51, 0x22, 0x26, 0xa1, 0xe5, 0x6f, 0xbe, 0x0b, 0x33, 0xaa, 0xb2, 0x08,
                    0x0a, 0xb4, 0x67, 0xd1, 0x4a, 0xad, 0xef, 0xf5, 0xdc, 0xd8, 0xad, 0xc6, 0x61,
                    0x3b, 0x92, 0x6b, 0xc9, 0x76, 0x01, 0xa4, 0xa1, 0xf1, 0x28, 0x77, 0x93, 0x75,
                    0x7b, 0x10, 0xd6, 0x8a, 0x93, 0x02, 0x30, 0x9a, 0x7a, 0x5c, 0x05, 0x62, 0x2b,
                    0xe1, 0xdb, 0xfb, 0xd4, 0xcd, 0x08, 0x50, 0x5f, 0x33, 0x37, 0xa4, 0x3e, 0x4f,
                    0xeb, 0x76, 0xac, 0x0f, 0xec, 0x48, 0xe0, 0x15, 0xd9, 0xc5, 0xa4, 0xdc, 0x1e,
                    0xed, 0xf0, 0x95, 0xfe, 0xe3, 0xbd, 0xc4, 0xb3, 0x73, 0x2c, 0x4c, 0xfd, 0x1b,
                    0x5b, 0x04, 0x9e, 0xc4, 0x34, 0xfe, 0xe2, 0x30, 0x23, 0x17, 0x94, 0xc2, 0xd0,
                    0xe6, 0x69, 0x74, 0x6d, 0x30, 0x36, 0x92, 0x64, 0xf8, 0xfc, 0x11, 0x29, 0x01,
                    0xcf, 0x11, 0xcb, 0xe1, 0xde, 0x80, 0x2f, 0x36, 0x9e, 0xf7, 0x61, 0xeb, 0xe1,
                    0xa0, 0x45, 0xca, 0xba, 0x36, 0x89, 0x1a, 0xab, 0x6e, 0x1c,
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
