//! ristretto255: a prime-order group built on edwards25519 (RFC 9496).
//!
//! edwards25519 is not of prime order: it has a cofactor of 8, which exposes a
//! small subgroup that complicates protocols. ristretto255 is a construction
//! that quotients that cofactor away, yielding a group of prime order `l`
//! (the [`super::Scalar`] field order) with:
//!
//! * a canonical 32-byte encoding ([`RistrettoPoint::compress`]) and decoding
//!   ([`RistrettoPoint::decompress`]) that rejects non-canonical inputs, so
//!   equal group elements always have equal encodings;
//! * equality testing that is independent of the internal representative
//!   ([`RistrettoPoint::ct_eq`] / `PartialEq`);
//! * the group operations (add, negate, scalar multiplication), which are the
//!   edwards25519 operations on a representative;
//! * a one-way map from 64 uniform bytes ([`RistrettoPoint::from_uniform_bytes`])
//!   suitable as the inner map of a hash-to-group construction.
//!
//! The implementation reuses the edwards25519 [`super::Point`] arithmetic and
//! the [`super::FieldElement`] `sqrt_ratio_m1` primitive; only the
//! encoding/decoding, equality, and Elligator map are specific to ristretto.

use super::{EdCurve, FieldElement, Point, Scalar};
use crate::curve::edwards::EdwardsCurve;
use crate::mp::ct::{Choice, CtEqual, CtOption, CtSelect, CtZero};
use std::ops::{Add, Mul, Neg, Sub};

// Field constants (big-endian), as specified by RFC 9496 Appendix A.
// `i = sqrt(-1)` is reused from the parent module as `FieldElement::SQRT_M1`.

/// `1 / sqrt(a - d)` with `a = -1`, i.e. `1 / sqrt(-1 - d)`.
const INVSQRT_A_MINUS_D: FieldElement = FieldElement::from_bytes_unchecked_be(&[
    0x78, 0x6c, 0x89, 0x05, 0xcf, 0xaf, 0xfc, 0xa2, 0x16, 0xc2, 0x7b, 0x91, 0xfe, 0x01, 0xd8, 0x40,
    0x9d, 0x2f, 0x16, 0x17, 0x5a, 0x41, 0x72, 0xbe, 0x99, 0xc8, 0xfd, 0xaa, 0x80, 0x5d, 0x40, 0xea,
]);

/// `sqrt(a*d - 1)` with `a = -1`, i.e. `sqrt(-d - 1)`.
const SQRT_AD_MINUS_ONE: FieldElement = FieldElement::from_bytes_unchecked_be(&[
    0x37, 0x69, 0x31, 0xbf, 0x2b, 0x83, 0x48, 0xac, 0x0f, 0x3c, 0xfc, 0xc9, 0x31, 0xf5, 0xd1, 0xfd,
    0xaf, 0x9d, 0x8e, 0x0c, 0x1b, 0x78, 0x54, 0xbd, 0x7e, 0x97, 0xf6, 0xa0, 0x49, 0x7b, 0x2e, 0x1b,
]);

/// `1 - d^2`.
const ONE_MINUS_D_SQ: FieldElement = FieldElement::from_bytes_unchecked_be(&[
    0x02, 0x90, 0x72, 0xa8, 0xb2, 0xb3, 0xe0, 0xd7, 0x99, 0x94, 0xab, 0xdd, 0xbe, 0x70, 0xdf, 0xe4,
    0x2c, 0x81, 0xa1, 0x38, 0xcd, 0x5e, 0x35, 0x0f, 0xe2, 0x7c, 0x09, 0xc1, 0x94, 0x5f, 0xc1, 0x76,
]);

/// `(d - 1)^2`.
const D_MINUS_ONE_SQ: FieldElement = FieldElement::from_bytes_unchecked_be(&[
    0x59, 0x68, 0xb3, 0x7a, 0xf6, 0x6c, 0x22, 0x41, 0x4c, 0xdc, 0xd3, 0x2f, 0x52, 0x9b, 0x4e, 0xeb,
    0xd2, 0x9e, 0x4a, 0x2c, 0xb0, 0x1e, 0x19, 0x99, 0x31, 0xad, 0x5a, 0xaa, 0x44, 0xed, 0x4d, 0x20,
]);

/// An element of the ristretto255 prime-order group.
///
/// Internally this wraps an edwards25519 [`Point`]; a given group element has
/// many valid representatives (differing by the 8-torsion), but they all
/// compress to the same 32-byte string and compare equal.
#[derive(Clone, Debug)]
pub struct RistrettoPoint(Point);

impl RistrettoPoint {
    /// The group identity element.
    pub const IDENTITY: Self = RistrettoPoint(Point::IDENTITY);

    /// The ristretto255 generator (the edwards25519 base point).
    pub const BASEPOINT: Self = RistrettoPoint(Point::GENERATOR);

    /// Encode this group element to its canonical 32-byte little-endian string.
    ///
    /// Equal group elements always produce identical bytes, regardless of the
    /// internal representative.
    pub fn compress(&self) -> [u8; 32] {
        let p = &self.0;
        let one = FieldElement::one();

        let u1 = &(&p.z + &p.y) * &(&p.z - &p.y);
        let u2 = &p.x * &p.y;
        let (_, invsqrt) = FieldElement::sqrt_ratio_m1(&one, &(&u1 * &u2.square()));
        let den1 = &invsqrt * &u1;
        let den2 = &invsqrt * &u2;
        let z_inv = &(&den1 * &den2) * &p.t;

        let ix = &p.x * &FieldElement::SQRT_M1;
        let iy = &p.y * &FieldElement::SQRT_M1;
        let enchanted_denominator = &den1 * &INVSQRT_A_MINUS_D;

        let rotate = (&p.t * &z_inv).is_negative_ct();
        let x = FieldElement::ct_select(rotate, &iy, &p.x);
        let y = FieldElement::ct_select(rotate, &ix, &p.y);
        let den_inv = FieldElement::ct_select(rotate, &enchanted_denominator, &den2);

        // conditionally negate y so that x is non-negative
        let y = FieldElement::ct_select((&x * &z_inv).is_negative_ct(), &(-&y), &y);

        let s = (&den_inv * &(&p.z - &y)).abs();
        s.to_bytes_le()
    }

    /// Decode a 32-byte string back into a group element.
    ///
    /// Returns a present [`CtOption`] only when `bytes` is a *canonical*
    /// ristretto255 encoding; non-canonical field encodings, negative `s`, and
    /// the other ill-formed cases of RFC 9496 are rejected.
    pub fn decompress(bytes: &[u8; 32]) -> CtOption<RistrettoPoint> {
        // s must be a canonical field element (< p) ...
        let s = match FieldElement::from_bytes_le(bytes) {
            Some(s) => s,
            None => return CtOption::from((Choice::FALSE, RistrettoPoint::IDENTITY)),
        };
        // ... and non-negative (even canonical encoding).
        if s.is_negative_ct().is_true() {
            return CtOption::from((Choice::FALSE, RistrettoPoint::IDENTITY));
        }

        let one = FieldElement::one();
        let ss = s.square();
        let u1 = &one - &ss; // 1 + a*s^2, a = -1
        let u2 = &one + &ss; // 1 - a*s^2
        let u2_sq = u2.square();

        // v = a*d*u1^2 - u2^2 = -d*u1^2 - u2^2
        let v = &(-&(&EdCurve::D * &u1.square())) - &u2_sq;

        let (was_square, invsqrt) = FieldElement::sqrt_ratio_m1(&one, &(&v * &u2_sq));
        let den_x = &invsqrt * &u2;
        let den_y = &(&invsqrt * &den_x) * &v;

        let x = (&(&s + &s) * &den_x).abs();
        let y = &u1 * &den_y;
        let t = &x * &y;

        let ok = was_square & t.is_negative_ct().negate() & y.ct_nonzero();
        CtOption::from((ok, RistrettoPoint(Point::from_affine(&x, &y))))
    }

    /// Constant-time equality, independent of the internal representative.
    pub fn ct_eq(&self, other: &RistrettoPoint) -> Choice {
        let p = &self.0;
        let q = &other.0;
        let x1y2 = &p.x * &q.y;
        let y1x2 = &p.y * &q.x;
        let x1x2 = &p.x * &q.x;
        let y1y2 = &p.y * &q.y;
        x1y2.ct_eq(&y1x2) | x1x2.ct_eq(&y1y2)
    }

    /// Group doubling (`self + self`).
    pub fn double(&self) -> RistrettoPoint {
        RistrettoPoint(self.0.double())
    }

    /// Scalar multiplication `k * self` (constant time).
    pub fn scale(&self, k: &Scalar) -> RistrettoPoint {
        RistrettoPoint(self.0.scale(k))
    }

    /// Fixed-base scalar multiplication `k * BASEPOINT` (constant time).
    pub fn mul_base(k: &Scalar) -> RistrettoPoint {
        RistrettoPoint(Point::mul_base(k))
    }

    /// The ristretto255 one-way map from 64 uniformly random bytes.
    ///
    /// Each 32-byte half is reduced to a field element (after masking the high
    /// bit, per RFC 9496) and run through the Elligator map; the two resulting
    /// points are added. With uniform input the output is statistically close to
    /// uniform over the group, making this the inner map of a hash-to-group
    /// construction (the outer domain-separated hash is left to the caller).
    pub fn from_uniform_bytes(bytes: &[u8; 64]) -> RistrettoPoint {
        let r0 = fe_from_wide_half(&bytes[0..32]);
        let r1 = fe_from_wide_half(&bytes[32..64]);
        RistrettoPoint(elligator(&r0).add(&elligator(&r1)))
    }
}

/// Reduce a 32-byte little-endian half (with its high bit masked) modulo p.
fn fe_from_wide_half(half: &[u8]) -> FieldElement {
    let mut wide = [0u8; FieldElement::SIZE_BYTES * 2];
    wide[..32].copy_from_slice(half);
    wide[31] &= 0x7f; // mask off bit 255, per RFC 9496
    FieldElement::init_from_wide_bytes_le(wide)
}

/// The ristretto255 Elligator map: a field element to an edwards25519 point.
fn elligator(t: &FieldElement) -> Point {
    let one = FieldElement::one();
    let minus_one = -&one;

    let r = &FieldElement::SQRT_M1 * &t.square();
    let u = &(&r + &one) * &ONE_MINUS_D_SQ;
    let v = &(&minus_one - &(&r * &EdCurve::D)) * &(&r + &EdCurve::D);

    let (was_square, s) = FieldElement::sqrt_ratio_m1(&u, &v);
    let s_prime = -&(&s * t).abs(); // -|s*t|
    let s = FieldElement::ct_select(was_square, &s, &s_prime);
    let n = FieldElement::ct_select(was_square, &minus_one, &r);

    let nn = &(&(&n * &(&r - &one)) * &D_MINUS_ONE_SQ) - &v;
    let ss = s.square();
    let w0 = &(&s + &s) * &v;
    let w1 = &nn * &SQRT_AD_MINUS_ONE;
    let w2 = &one - &ss;
    let w3 = &one + &ss;

    // extended coordinates (X:Y:Z:T)
    Point {
        x: &w0 * &w3,
        y: &w2 * &w1,
        z: &w1 * &w3,
        t: &w0 * &w2,
    }
}

impl PartialEq for RistrettoPoint {
    fn eq(&self, other: &Self) -> bool {
        self.ct_eq(other).is_true()
    }
}
impl Eq for RistrettoPoint {}

impl Neg for RistrettoPoint {
    type Output = RistrettoPoint;
    fn neg(self) -> RistrettoPoint {
        RistrettoPoint(-&self.0)
    }
}
impl<'a> Neg for &'a RistrettoPoint {
    type Output = RistrettoPoint;
    fn neg(self) -> RistrettoPoint {
        RistrettoPoint(-&self.0)
    }
}

impl<'a, 'b> Add<&'b RistrettoPoint> for &'a RistrettoPoint {
    type Output = RistrettoPoint;
    fn add(self, other: &'b RistrettoPoint) -> RistrettoPoint {
        RistrettoPoint(self.0.add(&other.0))
    }
}
impl<'a, 'b> Sub<&'b RistrettoPoint> for &'a RistrettoPoint {
    type Output = RistrettoPoint;
    fn sub(self, other: &'b RistrettoPoint) -> RistrettoPoint {
        RistrettoPoint(&self.0 - &other.0)
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a RistrettoPoint {
    type Output = RistrettoPoint;
    fn mul(self, k: &'b Scalar) -> RistrettoPoint {
        self.scale(k)
    }
}
impl<'a, 'b> Mul<&'b RistrettoPoint> for &'a Scalar {
    type Output = RistrettoPoint;
    fn mul(self, p: &'b RistrettoPoint) -> RistrettoPoint {
        p.scale(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hexval(c: u8) -> u8 {
        match c {
            b'0'..=b'9' => c - b'0',
            b'a'..=b'f' => c - b'a' + 10,
            _ => panic!("bad hex digit"),
        }
    }

    fn unhex<const N: usize>(s: &str) -> [u8; N] {
        let b = s.as_bytes();
        assert_eq!(b.len(), 2 * N, "hex string of wrong length");
        let mut out = [0u8; N];
        for (i, o) in out.iter_mut().enumerate() {
            *o = (hexval(b[2 * i]) << 4) | hexval(b[2 * i + 1]);
        }
        out
    }

    // RFC 9496 Appendix A.1: encodings of the first 16 multiples of the base point.
    const MULTIPLES: [&str; 16] = [
        "0000000000000000000000000000000000000000000000000000000000000000",
        "e2f2ae0a6abc4e71a884a961c500515f58e30b6aa582dd8db6a65945e08d2d76",
        "6a493210f7499cd17fecb510ae0cea23a110e8d5b901f8acadd3095c73a3b919",
        "94741f5d5d52755ece4f23f044ee27d5d1ea1e2bd196b462166b16152a9d0259",
        "da80862773358b466ffadfe0b3293ab3d9fd53c5ea6c955358f568322daf6a57",
        "e882b131016b52c1d3337080187cf768423efccbb517bb495ab812c4160ff44e",
        "f64746d3c92b13050ed8d80236a7f0007c3b3f962f5ba793d19a601ebb1df403",
        "44f53520926ec81fbd5a387845beb7df85a96a24ece18738bdcfa6a7822a176d",
        "903293d8f2287ebe10e2374dc1a53e0bc887e592699f02d077d5263cdd55601c",
        "02622ace8f7303a31cafc63f8fc48fdc16e1c8c8d234b2f0d6685282a9076031",
        "20706fd788b2720a1ed2a5dad4952b01f413bcf0e7564de8cdc816689e2db95f",
        "bce83f8ba5dd2fa572864c24ba1810f9522bc6004afe95877ac73241cafdab42",
        "e4549ee16b9aa03099ca208c67adafcafa4c3f3e4e5303de6026e3ca8ff84460",
        "aa52e000df2e16f55fb1032fc33bc42742dad6bd5a8fc0be0167436c5948501f",
        "46376b80f409b29dc2b5f6f0c52591990896e5716f41477cd30085ab7f10301e",
        "e0c418f7c8d9c4cdd7395b93ea124f3ad99021bb681dfc3302a9d99a2e53e64e",
    ];

    // RFC 9496 Appendix A.3: (64-byte uniform input, resulting encoding).
    const UNIFORM: [(&str, &str); 7] = [
        ("5d1be09e3d0c82fc538112490e35701979d99e06ca3e2b5b54bffe8b4dc772c14d98b696a1bbfb5ca32c436cc61c16563790306c79eaca7705668b47dffe5bb6",
         "3066f82a1a747d45120d1740f14358531a8f04bbffe6a819f86dfe50f44a0a46"),
        ("f116b34b8f17ceb56e8732a60d913dd10cce47a6d53bee9204be8b44f6678b270102a56902e2488c46120e9276cfe54638286b9e4b3cdb470b542d46c2068d38",
         "f26e5b6f7d362d2d2a94c5d0e7602cb4773c95a2e5c31a64f133189fa76ed61b"),
        ("8422e1bbdaab52938b81fd602effb6f89110e1e57208ad12d9ad767e2e25510c27140775f9337088b982d83d7fcf0b2fa1edffe51952cbe7365e95c86eaf325c",
         "006ccd2a9e6867e6a2c5cea83d3302cc9de128dd2a9a57dd8ee7b9d7ffe02826"),
        ("ac22415129b61427bf464e17baee8db65940c233b98afce8d17c57beeb7876c2150d15af1cb1fb824bbd14955f2b57d08d388aab431a391cfc33d5bafb5dbbaf",
         "f8f0c87cf237953c5890aec3998169005dae3eca1fbb04548c635953c817f92a"),
        ("165d697a1ef3d5cf3c38565beefcf88c0f282b8e7dbd28544c483432f1cec7675debea8ebb4e5fe7d6f6e5db15f15587ac4d4d4a1de7191e0c1ca6664abcc413",
         "ae81e7dedf20a497e10c304a765c1767a42d6e06029758d2d7e8ef7cc4c41179"),
        ("a836e6c9a9ca9f1e8d486273ad56a78c70cf18f0ce10abb1c7172ddd605d7fd2979854f47ae1ccf204a33102095b4200e5befc0465accc263175485f0e17ea5c",
         "e2705652ff9f5e44d3e841bf1c251cf7dddb77d140870d1ab2ed64f1a9ce8628"),
        ("2cdc11eaeb95daf01189417cdddbf95952993aa9cb9c640eb5058d09702c74622c9965a697a3b345ec24ee56335b556e677b30e6f90ac77d781064f866a3c982",
         "80bd07262511cdde4863f8a7434cef696750681cb9510eea557088f76d9e5065"),
    ];

    // RFC 9496 Appendix A.2: encodings that must be rejected, plus a
    // non-canonical (s >= p) field encoding (2^255 - 1).
    const BAD: [&str; 17] = [
        "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
        "0100000000000000000000000000000000000000000000000000000000000000",
        "ed57ffd8c914fb201471d1c3d245ce3c746fcbe63a3679d51b6a516ebebe0e20",
        "c34c4e1826e5d403b78e246e88aa051c36ccf0aafebffe137d148a2bf9104562",
        "c940e5a4404157cfb1628b108db051a8d439e1a421394ec4ebccb9ec92a8ac78",
        "f1c6165d33367351b0da8f6e4511010c68174a03b6581212c71c0e1d026c3c72",
        "87260f7a2f12495118360f02c26a470f450dadf34a413d21042b43b9d93e1309",
        "26948d35ca62e643e26a831777332e6bafeb9d08e4268b650f1f5bbd8d81d371",
        "de6a7b00deadc788eb6b6c8d20c0ae96c2f2019078fa604fee5b87d6e989ad7b",
        "f4a9e534fc0d216c44b218fa0c42d99635a0127ee2e53c712f70609649fdff22",
        "2810e5cbc2cc4d4eece54f61c6f69758e289aa7ab440b3cbeaa21995c2f4232b",
        "3eb858e78f5a7254d8c9731174a94f76755fd3941c0ac93735c07ba14579630e",
        "a45fdc55c76448c049a1ab33f17023edfb2be3581e9c7aade8a6125215e04220",
        "d483fe813c6ba647ebbfd3ec41adca1c6130c2beeee9d9bf065c8d151c5f396e",
        "8a2e1d30050198c65a54483123960ccc38aef6848e1ec8f5f780e8523769ba32",
        "32888462f8b486c68ad7dd9610be5192bbeaf3b443951ac1a8118419d9fa097b",
        "5c37cc491da847cfeb9281d407efc41e15144c876e0170b499a96a22ed31e01e",
    ];

    #[test]
    fn basepoint_multiples_kat() {
        let mut acc = RistrettoPoint::IDENTITY;
        for (k, expected) in MULTIPLES.iter().enumerate() {
            assert_eq!(acc.compress(), unhex::<32>(expected), "[{}]B encoding", k);
            acc = &acc + &RistrettoPoint::BASEPOINT;
        }
    }

    #[test]
    fn decompress_roundtrip() {
        let mut acc = RistrettoPoint::IDENTITY;
        for (k, enc) in MULTIPLES.iter().enumerate() {
            let bytes = unhex::<32>(enc);
            let p = RistrettoPoint::decompress(&bytes)
                .into_option()
                .unwrap_or_else(|| panic!("[{}]B should decompress", k));
            // decoded value matches the independently accumulated point ...
            assert_eq!(p, acc, "[{}]B decode mismatch", k);
            // ... and re-encodes to the same bytes.
            assert_eq!(p.compress(), bytes, "[{}]B reencode", k);
            acc = &acc + &RistrettoPoint::BASEPOINT;
        }
    }

    #[test]
    fn bad_encodings_rejected() {
        for enc in BAD.iter() {
            let bytes = unhex::<32>(enc);
            assert!(
                RistrettoPoint::decompress(&bytes).into_option().is_none(),
                "must reject {}",
                enc
            );
        }
    }

    #[test]
    fn uniform_bytes_kat() {
        for (input, expected) in UNIFORM.iter() {
            let p = RistrettoPoint::from_uniform_bytes(&unhex::<64>(input));
            assert_eq!(p.compress(), unhex::<32>(expected), "input {}", input);
        }
    }

    #[test]
    fn group_law() {
        let b = RistrettoPoint::BASEPOINT;
        // identity behaviour
        assert_eq!(&b + &RistrettoPoint::IDENTITY, b);
        assert_eq!(&b - &b, RistrettoPoint::IDENTITY);
        // doubling and small scalars agree with repeated addition
        assert_eq!(b.double(), &b + &b);
        let two = Scalar::from_u64(2);
        let three = Scalar::from_u64(3);
        assert_eq!(b.scale(&two), b.double());
        assert_eq!(b.scale(&three), &b.double() + &b);
        // fixed-base agrees with variable-base
        assert_eq!(RistrettoPoint::mul_base(&three), b.scale(&three));
        // negation
        assert_eq!(&b + &(-&b), RistrettoPoint::IDENTITY);
        // a scalar multiple is not the identity (sanity)
        assert_ne!(b.scale(&three), RistrettoPoint::IDENTITY);
    }

    #[test]
    fn equality_independent_of_representative() {
        // The accumulated [3]B has a non-trivial Z (from repeated addition),
        // while the decompressed one is in (x, y, 1, xy) form. They must still
        // compare equal.
        let three = Scalar::from_u64(3);
        let accumulated = RistrettoPoint::BASEPOINT.scale(&three);
        let decoded = RistrettoPoint::decompress(&unhex::<32>(MULTIPLES[3]))
            .into_option()
            .unwrap();
        assert_eq!(accumulated, decoded);
        assert!(accumulated.ct_eq(&decoded).is_true());
    }
}
