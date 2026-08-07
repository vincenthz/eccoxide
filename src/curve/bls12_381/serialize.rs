//! Standard BLS12-381 point encodings
//!
//! The encoding introduced by zcash and reused by the IETF BLS signature
//! specification (see
//! <https://github.com/zkcrypto/bls12_381/blob/main/src/notes/serialization.rs>).
//!
//! All serialization are big-endian and fixed size.
//!
//! It currently exists in two flavours:
//!
//! * *uncompressed*: `x || y`, two coordinates,
//! * *compressed*: `x` only, with one bit telling which of the two square roots
//!   `y` is.
//!
//! A coordinate of `G1` is an [`Fp`] element, encoded as 48 big-endian bytes; a
//! coordinate of `G2` is an [`Fp2`] element `c0 + c1*u`, encoded as the 96 bytes
//! `c1 || c0` — the *imaginary* part first. So the four sizes are 48/96 bytes
//! for `G1` and 96/192 bytes for `G2`.
//!
//! Since `p < 2^381`, the top three bits of the leading byte of a canonical
//! coordinate are always clear, and the format stores its three flags there:
//!
//! | bit | flag | meaning |
//! |-----|------|---------|
//! | 7 | [`COMPRESSION_FLAG`] | the encoding carries only `x` |
//! | 6 | [`INFINITY_FLAG`] | the point is the identity; all other bits are zero |
//! | 5 | [`SORT_FLAG`] | (compressed only) `y` is the larger of the two roots |
//!
//! "Larger" is the lexicographic order on the canonical representative: `y` is
//! the larger root when `y > (p - 1)/2`, and for an `Fp2` coordinate the `c1`
//! component is compared first (matching the `c1 || c0` byte order). Note this
//! is *not* the parity convention of [`Sign`](crate::curve::field::Sign), which
//! the generic [`compress`](crate::curve::affine::Point::compress) /
//! [`decompress`](crate::curve::affine::Point::decompress) pair uses, hence the
//! sign fixup when reading a compressed encoding.
//!
//! The identity has a single encoding per flavour: the infinity flag set (plus
//! the compression flag when compressed) and every other bit zero.
//!
//! The coordinates themselves are written by [`Fp::to_bytes`] and
//! [`Fp2::to_bytes`], which are in this very order; what this module adds is
//! the flag bits and the ordering the sort flag is defined by, and
//! [`bls12_381_define_point_serialization`] assembles them into the methods of
//! a group's point types.

use super::fp::Fp;
use super::fp2::Fp2;
use crate::mp::ct::{Choice, CtLesser, CtZero};
use crate::params::bls12_381::P_MINUS1_DIV2_BYTES;

/// Bit 7 of the leading byte: only the x-coordinate is present.
pub(crate) const COMPRESSION_FLAG: u8 = 0b1000_0000;

/// Bit 6 of the leading byte: the encoded point is the identity.
pub(crate) const INFINITY_FLAG: u8 = 0b0100_0000;

/// Bit 5 of the leading byte: in a compressed encoding, the dropped `y` is the
/// lexicographically larger of the two square roots.
pub(crate) const SORT_FLAG: u8 = 0b0010_0000;

/// The three flag bits of the leading byte.
pub(super) const FLAGS_MASK: u8 = COMPRESSION_FLAG | INFINITY_FLAG | SORT_FLAG;

/// The validated flags of a compressed encoding.
pub(crate) enum Compressed {
    Infinity,
    /// A point whose dropped `y` is the larger root when `sort` is true.
    Point {
        sort: Choice,
    },
}

/// Fold `flag` into `byte` when `set`, without branching on the choice.
fn set_flag(byte: &mut u8, flag: u8, set: Choice) {
    *byte |= flag & 0u8.wrapping_sub(set.to_u1());
}

/// Whether every bit of the encoding outside the flags is zero, as the identity
/// encoding requires.
fn is_zero_payload(bytes: &[u8]) -> bool {
    let (first, rest) = bytes.split_first().expect("non-empty encoding");
    rest.iter().fold(first & !FLAGS_MASK, |acc, b| acc | b) == 0
}

/// The encoding of the identity: a leading byte of flags and nothing else.
const fn infinity<const N: usize>(flags: u8) -> [u8; N] {
    let mut out = [0u8; N];
    out[0] = flags;
    out
}

/// The compressed encoding of the identity.
pub(crate) const fn infinity_compressed<const N: usize>() -> [u8; N] {
    infinity(COMPRESSION_FLAG | INFINITY_FLAG)
}

/// The uncompressed encoding of the identity.
pub(crate) const fn infinity_uncompressed<const N: usize>() -> [u8; N] {
    infinity(INFINITY_FLAG)
}

/// Set the flags of a compressed encoding whose x-coordinate is already written.
pub(crate) fn write_compressed_flags(out: &mut [u8], y_is_largest: Choice) {
    out[0] |= COMPRESSION_FLAG;
    set_flag(&mut out[0], SORT_FLAG, y_is_largest);
}

/// Validate the flags of a compressed encoding.
pub(crate) fn read_compressed_flags(bytes: &[u8]) -> Option<Compressed> {
    let flags = bytes[0] & FLAGS_MASK;

    if flags & COMPRESSION_FLAG == 0 {
        return None;
    }
    if flags & INFINITY_FLAG != 0 {
        // no root to sort, so the sort flag must be clear as well
        if flags & SORT_FLAG != 0 || !is_zero_payload(bytes) {
            return None;
        }
        return Some(Compressed::Infinity);
    }
    Some(Compressed::Point {
        sort: Choice::from(flags & SORT_FLAG != 0),
    })
}

/// Validate the flags of an uncompressed encoding, and report whether they mark
/// the identity.
pub(crate) fn read_uncompressed_flags(bytes: &[u8]) -> Option<bool> {
    let flags = bytes[0] & FLAGS_MASK;

    // this flavour carries both coordinates, and the sort flag has no meaning
    if flags & (COMPRESSION_FLAG | SORT_FLAG) != 0 {
        return None;
    }
    if flags & INFINITY_FLAG != 0 {
        return is_zero_payload(bytes).then_some(true);
    }
    Some(false)
}

impl Fp {
    /// Whether this is the larger of `self` and `-self`, i.e. `self > (p - 1)/2`.
    #[inline(always)]
    pub(crate) fn is_largest(&self) -> Choice {
        let bytes = self.to_bytes();
        <&[u8; Fp::SIZE_BYTES]>::ct_lt(&P_MINUS1_DIV2_BYTES, &bytes)
    }
}

impl Fp2 {
    /// Whether this is the larger of `self` and `-self`, comparing `c1` first to
    /// match the `c1 || c0` byte order of [`Fp2::to_bytes`].
    #[inline(always)]
    pub(crate) fn is_largest(&self) -> Choice {
        self.c1.is_largest() | (self.c1.ct_zero() & self.c0.is_largest())
    }
}

/// Emit the standard encoding methods for one of the BLS12-381 groups.
///
/// It extends:
///
/// * `PointAffine`
/// * `Point`
///
/// with from_uncompressed / to_uncompressed / from_compressed / to_compressed
#[doc(hidden)]
#[macro_export]
macro_rules! bls12_381_define_point_serialization {
    ($FE:ident) => {
        /// Decode the x-coordinate of a compressed encoding whose flags are
        /// already validated, and recover the root `sort` asks for.
        ///
        /// Constant time: a rejection, whether of a non-canonical coordinate or
        /// of an x that is not on the curve, is carried as the choice of the
        /// returned [`CtOption`](crate::mp::ct::CtOption) rather than taken as a
        /// branch, and the carried point is then a placeholder.
        fn read_compressed_affine(
            bytes: &[u8; PointAffine::COMPRESSED_SIZE],
            sort: $crate::mp::ct::Choice,
        ) -> $crate::mp::ct::CtOption<PointAffine> {
            use $crate::curve::bls12_381::serialize;
            use $crate::mp::ct::{CtEqual, CtOption, CtSelect};

            // the flags share the leading byte with the x-coordinate; clear them
            // to recover its encoding, which then has to be canonical on its own
            let mut buf = *bytes;
            buf[0] &= !serialize::FLAGS_MASK;
            let (canonical, x) = $FE::from_bytes_ct(&buf).into_parts();

            // decompression picks between the two roots by parity, so swap to
            // the one the sort flag asks for
            let (on_curve, p) = affine::Point::decompress::<Curve>(&x, Sign::Positive).into_parts();
            let neg_y = -p.y.clone();
            let y = $FE::ct_select(p.y.is_largest().ct_eq(&sort), &p.y, &neg_y);

            CtOption::from((canonical & on_curve, PointAffine(affine::Point { x, y })))
        }

        /// Decode the coordinates of an uncompressed encoding whose flags are
        /// already validated.
        ///
        /// Constant time in the same way as [`read_compressed_affine`].
        fn read_uncompressed_affine(
            bytes: &[u8; PointAffine::UNCOMPRESSED_SIZE],
        ) -> $crate::mp::ct::CtOption<PointAffine> {
            use core::convert::TryInto;
            use $crate::mp::ct::CtOption;

            // with all three flags clear, the coordinates are the bytes as they
            // stand; the halves are of the right size by construction
            let (xb, yb) = bytes.split_at($FE::SIZE_BYTES);
            let (x_canonical, x) = $FE::from_bytes_ct(xb.try_into().unwrap()).into_parts();
            let (y_canonical, y) = $FE::from_bytes_ct(yb.try_into().unwrap()).into_parts();

            let (on_curve, p) = affine::Point::from_coordinate_ct::<Curve>(&x, &y).into_parts();
            CtOption::from((x_canonical & y_canonical & on_curve, PointAffine(p)))
        }

        impl PointAffine {
            /// Size in bytes of the compressed encoding.
            pub const COMPRESSED_SIZE: usize = $FE::SIZE_BYTES;

            /// Size in bytes of the uncompressed encoding.
            pub const UNCOMPRESSED_SIZE: usize = 2 * $FE::SIZE_BYTES;

            /// Serialize in the standard compressed format: the big-endian
            /// x-coordinate, with the compression flag and the bit telling
            /// which of the two roots `y` is in the top bits of the leading
            /// byte.
            ///
            /// See [`Self::from_compressed`] for the details of the format.
            pub fn to_compressed(&self) -> [u8; Self::COMPRESSED_SIZE] {
                use $crate::curve::bls12_381::serialize;

                let (x, y) = self.to_coordinate();
                let mut out = [0u8; Self::COMPRESSED_SIZE];
                x.to_slice(&mut out);
                serialize::write_compressed_flags(&mut out, y.is_largest());
                out
            }

            /// Serialize in the standard uncompressed format: the big-endian
            /// `x || y` coordinates.
            ///
            /// All flag bits are clear (no compression/infinity/sort).
            ///
            /// See [`Self::from_uncompressed`] for the details of the format.
            pub fn to_uncompressed(&self) -> [u8; Self::UNCOMPRESSED_SIZE] {
                let (x, y) = self.to_coordinate();
                let mut out = [0u8; Self::UNCOMPRESSED_SIZE];
                let (xb, yb) = out.split_at_mut($FE::SIZE_BYTES);
                x.to_slice(xb);
                y.to_slice(yb);
                out
            }

            /// Deserialize the standard compressed format, as produced by
            /// [`Self::to_compressed`].
            ///
            /// The encoding is the x-coordinate in big-endian form (for a `G2`
            /// coordinate `c0 + c1*u`, as `c1 || c0`) with three flags in the
            /// top bits of the leading byte, which a canonical coordinate
            /// leaves free: bit 7 marks the encoding compressed and must be
            /// set here, bit 6 marks the identity, and bit 5 tells whether the
            /// dropped `y` is the lexicographically larger of the two roots.
            ///
            /// `None` is returned unless the encoding is the canonical one of
            /// a point of the curve: the x-coordinate must be the canonical
            /// representative of its class (less than `p`), `x^3 + b` must be
            /// a square, and the flags must be consistent — in particular the
            /// identity, which this type cannot represent, is rejected rather
            /// than reported. Beware that lying on the curve does not imply
            /// membership of the prime-order subgroup, which is not checked.
            pub fn from_compressed(bytes: &[u8; Self::COMPRESSED_SIZE]) -> Option<Self> {
                use $crate::curve::bls12_381::serialize::{read_compressed_flags, Compressed};

                match read_compressed_flags(bytes)? {
                    // the identity has no affine representation
                    Compressed::Infinity => None,
                    Compressed::Point { sort } => read_compressed_affine(bytes, sort).into_option(),
                }
            }

            /// Deserialize the standard uncompressed format, as produced by
            /// [`Self::to_uncompressed`].
            ///
            /// The encoding is `x || y`, each coordinate in big-endian form
            /// (for a `G2` coordinate `c0 + c1*u`, as `c1 || c0`), with the
            /// flags of [`Self::from_compressed`] in the top bits of the
            /// leading byte: the compression and sort bits must be clear here,
            /// and bit 6 marks the identity.
            ///
            /// `None` is returned unless the encoding is the canonical one of
            /// a point of the curve: both coordinates must be canonical, they
            /// must satisfy the curve equation, and the flags must be
            /// consistent — in particular the identity, which this type cannot
            /// represent, is rejected rather than reported. Beware that lying
            /// on the curve does not imply membership of the prime-order
            /// subgroup, which is not checked.
            pub fn from_uncompressed(bytes: &[u8; Self::UNCOMPRESSED_SIZE]) -> Option<Self> {
                use $crate::curve::bls12_381::serialize::read_uncompressed_flags;

                if read_uncompressed_flags(bytes)? {
                    // the identity has no affine representation
                    None
                } else {
                    read_uncompressed_affine(bytes).into_option()
                }
            }
        }

        impl Point {
            /// Size in bytes of the compressed encoding.
            pub const COMPRESSED_SIZE: usize = PointAffine::COMPRESSED_SIZE;

            /// Size in bytes of the uncompressed encoding.
            pub const UNCOMPRESSED_SIZE: usize = PointAffine::UNCOMPRESSED_SIZE;

            /// Serialize in the standard compressed format, the identity
            /// included.
            ///
            /// See [`PointAffine::from_compressed`] for the details of the
            /// format. Note this normalizes the point, so it costs a field
            /// inversion.
            pub fn to_compressed(&self) -> [u8; Self::COMPRESSED_SIZE] {
                use $crate::curve::bls12_381::serialize;
                use $crate::mp::ct::CtSelect;

                // the identity has no affine coordinates, but the constant-time
                // conversion still hands back a placeholder to encode; the
                // result is then swapped for the identity encoding through a
                // branch-free select rather than by testing the point
                let (not_infinity, p) = self.to_affine_ct().into_parts();
                let infinity = serialize::infinity_compressed();
                CtSelect::ct_select(not_infinity, &p.to_compressed(), &infinity)
            }

            /// Serialize in the standard uncompressed format, the identity
            /// included.
            ///
            /// See [`PointAffine::from_uncompressed`] for the details of the
            /// format. Note this normalizes the point, so it costs a field
            /// inversion.
            pub fn to_uncompressed(&self) -> [u8; Self::UNCOMPRESSED_SIZE] {
                use $crate::curve::bls12_381::serialize;
                use $crate::mp::ct::CtSelect;

                // see `to_compressed` about the placeholder and the select
                let (not_infinity, p) = self.to_affine_ct().into_parts();
                let infinity = serialize::infinity_uncompressed();
                CtSelect::ct_select(not_infinity, &p.to_uncompressed(), &infinity)
            }

            /// Deserialize the standard compressed format, as produced by
            /// [`Self::to_compressed`].
            ///
            /// Same as [`PointAffine::from_compressed`], which documents the
            /// format, except that the encoding of the identity is accepted
            /// and gives [`Self::INFINITY`].
            pub fn from_compressed(bytes: &[u8; Self::COMPRESSED_SIZE]) -> Option<Self> {
                use $crate::curve::bls12_381::serialize::{read_compressed_flags, Compressed};

                match read_compressed_flags(bytes)? {
                    Compressed::Infinity => Some(Point::INFINITY),
                    Compressed::Point { sort } => read_compressed_affine(bytes, sort)
                        .map(Point::from)
                        .into_option(),
                }
            }

            /// Deserialize the standard uncompressed format, as produced by
            /// [`Self::to_uncompressed`].
            ///
            /// Same as [`PointAffine::from_uncompressed`], which documents the
            /// format, except that the encoding of the identity is accepted
            /// and gives [`Self::INFINITY`].
            pub fn from_uncompressed(bytes: &[u8; Self::UNCOMPRESSED_SIZE]) -> Option<Self> {
                use $crate::curve::bls12_381::serialize::read_uncompressed_flags;

                if read_uncompressed_flags(bytes)? {
                    Some(Point::INFINITY)
                } else {
                    read_uncompressed_affine(bytes)
                        .map(Point::from)
                        .into_option()
                }
            }
        }
    };
}

/// Emit the group-agnostic [standard encoding](self) unit tests.
///
/// Invoke from inside the `tests` module of a group; the group-specific
/// known-answer vectors stay in that module.
#[doc(hidden)]
#[macro_export]
macro_rules! bls12_381_point_serialization_unittest {
    () => {
        mod serialization {
            use super::super::{Point, PointAffine, Scalar};
            use $crate::params::bls12_381::P_BYTES;

            const COMPRESSION_FLAG: u8 = 0b1000_0000;
            const INFINITY_FLAG: u8 = 0b0100_0000;
            const SORT_FLAG: u8 = 0b0010_0000;

            /// A handful of points of the group to run the round-trips over.
            fn samples() -> Vec<PointAffine> {
                [1u64, 2, 3, 12345, 0xdead_beef, 0x1234_5678_9abc_def0]
                    .iter()
                    .map(|v| {
                        (Point::GENERATOR * &Scalar::from_u64(*v))
                            .to_affine()
                            .expect("a multiple of the generator is not the identity")
                    })
                    .collect()
            }

            /// Both flavours round-trip, and the affine and projective sides
            /// agree on the bytes.
            #[test]
            fn roundtrip() {
                for p in samples() {
                    let q = Point::from(&p);

                    let c = p.to_compressed();
                    assert_eq!(PointAffine::from_compressed(&c).unwrap(), p);
                    assert_eq!(q.to_compressed(), c);
                    assert_eq!(Point::from_compressed(&c).unwrap(), q);

                    let u = p.to_uncompressed();
                    assert_eq!(PointAffine::from_uncompressed(&u).unwrap(), p);
                    assert_eq!(q.to_uncompressed(), u);
                    assert_eq!(Point::from_uncompressed(&u).unwrap(), q);
                }
            }

            /// The compressed encoding is the x-coordinate of the uncompressed
            /// one, flag bits aside.
            #[test]
            fn compressed_is_the_uncompressed_x() {
                for p in samples() {
                    let mut c = p.to_compressed();
                    c[0] &= !(COMPRESSION_FLAG | INFINITY_FLAG | SORT_FLAG);
                    assert_eq!(c[..], p.to_uncompressed()[..PointAffine::COMPRESSED_SIZE]);
                }
            }

            /// Negating a point keeps its x-coordinate and flips the sort flag,
            /// since `y` and `-y` are on either side of `(p - 1)/2`.
            #[test]
            fn negation_flips_the_sort_flag() {
                for p in samples() {
                    let neg = (-&Point::from(&p)).to_affine().unwrap();
                    let (a, b) = (p.to_compressed(), neg.to_compressed());
                    assert_eq!(a[0] ^ b[0], SORT_FLAG);
                    assert_eq!(a[1..], b[1..]);
                }
            }

            /// The identity has one encoding per flavour, which only the
            /// projective type can represent.
            #[test]
            fn infinity_roundtrip() {
                let inf = Point::INFINITY;

                let c = inf.to_compressed();
                assert_eq!(c[0], COMPRESSION_FLAG | INFINITY_FLAG);
                assert!(c[1..].iter().all(|b| *b == 0));
                assert_eq!(Point::from_compressed(&c).unwrap(), inf);
                assert!(PointAffine::from_compressed(&c).is_none());

                let u = inf.to_uncompressed();
                assert_eq!(u[0], INFINITY_FLAG);
                assert!(u[1..].iter().all(|b| *b == 0));
                assert_eq!(Point::from_uncompressed(&u).unwrap(), inf);
                assert!(PointAffine::from_uncompressed(&u).is_none());
            }

            /// Any projective representation of the identity encodes to the
            /// identity, not just the canonical one: the encoders never look at
            /// the placeholder coordinates the constant-time affine conversion
            /// hands back for it.
            #[test]
            fn infinity_encoding_ignores_the_representation() {
                let g = Point::GENERATOR;
                // `(0 : y : 0)` with a `y` of its own, rather than `Point::INFINITY`
                let inf = &g - &g;
                assert_eq!(inf, Point::INFINITY);
                assert_eq!(inf.to_compressed(), Point::INFINITY.to_compressed());
                assert_eq!(inf.to_uncompressed(), Point::INFINITY.to_uncompressed());
            }

            /// Every inconsistent flag combination is rejected.
            #[test]
            fn rejects_bad_flags() {
                let p = PointAffine::GENERATOR;
                let (c, u) = (p.to_compressed(), p.to_uncompressed());

                // a compressed encoding must say so, an uncompressed one must not
                let mut b = c;
                b[0] &= !COMPRESSION_FLAG;
                assert!(Point::from_compressed(&b).is_none());
                let mut b = u;
                b[0] |= COMPRESSION_FLAG;
                assert!(Point::from_uncompressed(&b).is_none());

                // the sort flag is only part of the compressed flavour
                let mut b = u;
                b[0] |= SORT_FLAG;
                assert!(Point::from_uncompressed(&b).is_none());

                // the identity flag with a non-zero payload, in both flavours
                let mut b = c;
                b[0] |= INFINITY_FLAG;
                assert!(Point::from_compressed(&b).is_none());
                let mut b = u;
                b[0] |= INFINITY_FLAG;
                assert!(Point::from_uncompressed(&b).is_none());

                // an identity encoding with a stray bit in the payload
                let mut b = Point::INFINITY.to_compressed();
                b[Point::COMPRESSED_SIZE - 1] = 1;
                assert!(Point::from_compressed(&b).is_none());
                let mut b = Point::INFINITY.to_uncompressed();
                b[Point::UNCOMPRESSED_SIZE - 1] = 1;
                assert!(Point::from_uncompressed(&b).is_none());

                // the identity has no root to sort
                let mut b = Point::INFINITY.to_compressed();
                b[0] |= SORT_FLAG;
                assert!(Point::from_compressed(&b).is_none());
            }

            /// A coordinate that is not the canonical representative of its
            /// class (here `p` itself) is rejected, in either position.
            #[test]
            fn rejects_non_canonical_coordinate() {
                // the leading component of x, in both groups
                let mut b = [0u8; PointAffine::COMPRESSED_SIZE];
                b[..P_BYTES.len()].copy_from_slice(&P_BYTES);
                b[0] |= COMPRESSION_FLAG;
                assert!(Point::from_compressed(&b).is_none());

                // and the leading component of y
                let half = PointAffine::UNCOMPRESSED_SIZE / 2;
                let mut b = PointAffine::GENERATOR.to_uncompressed();
                b[half..half + P_BYTES.len()].copy_from_slice(&P_BYTES);
                assert!(Point::from_uncompressed(&b).is_none());
            }

            /// An x-coordinate for which `x^3 + b` is not a square has no
            /// point to decompress to. Perturbing a valid x hits such an x
            /// about half the time, so a handful of tries is plenty; whatever
            /// is accepted has to round-trip.
            #[test]
            fn rejects_x_off_the_curve() {
                let base = PointAffine::GENERATOR.to_compressed();
                let mut rejected = 0;
                for i in 1..=32u8 {
                    let mut b = base;
                    b[PointAffine::COMPRESSED_SIZE - 1] ^= i;
                    match PointAffine::from_compressed(&b) {
                        None => rejected += 1,
                        Some(q) => assert_eq!(q.to_compressed(), b),
                    }
                }
                assert!(rejected > 0, "no perturbation of x fell off the curve");
            }

            /// Likewise, a `y` that does not match the x-coordinate is rejected
            /// by the uncompressed flavour.
            #[test]
            fn rejects_y_off_the_curve() {
                let mut b = PointAffine::GENERATOR.to_uncompressed();
                b[PointAffine::UNCOMPRESSED_SIZE - 1] ^= 1;
                assert!(Point::from_uncompressed(&b).is_none());
            }
        }
    };
}

#[cfg(test)]
mod tests {
    use crate::curve::bls12_381::{Fp, Fp2};
    use crate::curve::field::Field;

    /// The `is_largest` predicates split each pair `{y, -y}` in two: exactly one
    /// of the two is the larger, which is what makes one bit enough to tell them
    /// apart. Zero is its own negation and is not the larger of anything.
    #[test]
    fn exactly_one_root_is_the_largest() {
        for v in [1u64, 2, 3, 4, 5, 0xdead_beef, 0xffff_ffff_ffff_ffff] {
            let a = Fp::from_u64(v);
            assert_ne!(
                a.is_largest().is_true(),
                (-&a).is_largest().is_true(),
                "both or neither of +/-{} are the largest",
                v
            );
        }
        assert!(!Fp::ZERO.is_largest().is_true());
        assert!(!Fp2::ZERO.is_largest().is_true());

        // an Fp2 element has three shapes to cover: both components non-zero,
        // and either one of them zero, the c1 == 0 case being the one where the
        // comparison falls through to c0
        for (c0, c1) in [(2u64, 3), (7, 0), (0, 11), (1, 1)] {
            let a = Fp2::new(Fp::from_u64(c0), Fp::from_u64(c1));
            assert_ne!(
                a.is_largest().is_true(),
                (-&a).is_largest().is_true(),
                "both or neither of +/-({} + {}u) are the largest",
                c0,
                c1
            );
        }
    }

    /// Small values sit below `(p - 1)/2` and their negations above, and an
    /// `Fp2` element whose `c1` is zero is ordered by its `c0` alone.
    #[test]
    fn largest_is_the_half_above_the_midpoint() {
        let one = Fp::from_u64(1);
        assert!(!one.is_largest().is_true());
        assert!((-&one).is_largest().is_true()); // p - 1

        for v in [1u64, 42, 0xdead_beef] {
            let c0 = Fp::from_u64(v);
            let real = Fp2::new(c0.clone(), Fp::ZERO);
            assert_eq!(real.is_largest().is_true(), c0.is_largest().is_true());
            assert_eq!(
                (-&real).is_largest().is_true(),
                (-&c0).is_largest().is_true()
            );
        }
    }

    /// The `c1` component decides the order regardless of `c0`.
    #[test]
    fn largest_compares_c1_first() {
        let small = Fp::from_u64(1);
        let large = -Fp::from_u64(1); // p - 1, above the midpoint
        assert!(!Fp2::new(large.clone(), small.clone())
            .is_largest()
            .is_true());
        assert!(Fp2::new(small, large).is_largest().is_true());
    }
}
