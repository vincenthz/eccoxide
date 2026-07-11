//! ECDSA signature scheme (SEC1 / FIPS 186-4) over the SEC2 weierstrass curves.
//!
//! [`sign`] and [`verify`] are generic over a single [`EcdsaOperations`]
//! trait bundling everything a concrete instantiation needs:
//!
//! * the curve group it operates on, as the associated
//!   [`Point`](EcdsaOperations::Point) type — a [`CurveGroup`] (with the
//!   [`Field`] of its scalars) providing the curve arithmetic;
//! * [`x_mod_n`](EcdsaOperations::x_mod_n), the affine x-coordinate of a point
//!   reduced modulo the group order, the one extra curve operation ECDSA needs
//!   beyond the generic group law;
//! * [`hash_to_scalar`](EcdsaOperations::hash_to_scalar), the message-to-scalar
//!   hashing: digest a message with a specific hash function (SHA-256, ...)
//!   and convert the digest to a scalar of the curve.
//!
//! Marker types implementing it are provided pairing every enabled curve with
//! the SHA-2 family (e.g. [`P256R1_Sha256`] for ECDSA on NIST P-256 with
//! SHA-256).
//!
//! A signature is the generic [`Signature`] pair of scalars; fixed size
//! serialization is provided per curve (`to_bytes` / `from_bytes`).
//!
//! ```
//! # #[cfg(all(feature = "ecdsa", feature = "p256r1"))] {
//! use eccoxide::protocol::ecdsa::{self, P256R1_Sha256};
//! use eccoxide::curve::sec2::p256r1::{Point, Scalar};
//!
//! # let (secret_bytes, nonce_bytes) = ([3u8; 64], [7u8; 64]);
//! // secret key, and a fresh random nonce for this one signature
//! let secret = Scalar::init_from_wide_bytes_be(secret_bytes);
//! let nonce = Scalar::init_from_wide_bytes_be(nonce_bytes);
//!
//! let public: Point = ecdsa::public_key(&secret);
//! let signature = ecdsa::sign::<P256R1_Sha256>(&secret, &nonce, b"message")
//!     .into_option()
//!     .unwrap();
//! assert!(ecdsa::verify::<P256R1_Sha256>(&public, b"message", &signature));
//! # }
//! ```
//!
//! ## Custom hash
//!
//! To use a hash function not provided here, implement [`EcdsaOperations`] on
//! a marker type: set its [`Point`](EcdsaOperations::Point) (and matching
//! [`Scalar`](EcdsaOperations::Scalar)) to the curve's, delegate
//! [`x_mod_n`](EcdsaOperations::x_mod_n) to the curve's own (e.g.
//! [`p256r1::x_mod_n`]), and implement
//! [`hash_to_scalar`](EcdsaOperations::hash_to_scalar) by digesting the
//! message and converting the digest with the curve's `bits2int` transform
//! (e.g. [`p256r1::digest_to_scalar`]), which keeps the leftmost
//! `min(qlen, 8*N)` bits of the digest and reduces them modulo the curve
//! order — so any digest size can be paired with any curve. Signing a digest
//! computed externally is also possible directly with [`sign_hashed`] and
//! [`verify_hashed`] after the same conversion.
//!
//! ## Nonce
//!
//! ECDSA requires a secret nonce `k` that is unique and uniformly random for
//! every signature. Reusing a nonce for two different messages, or using a
//! biased nonce, leaks the secret key. This module takes the nonce as a
//! parameter and does not generate it: use either fresh randomness from a
//! CSPRNG (e.g. `Scalar::init_from_wide_bytes_be` over `2 * SCALAR_SIZE`
//! random bytes, to avoid modulo bias) or a deterministic derivation such as
//! RFC 6979.

use crate::curve::field::Field;
use crate::curve::CurveGroup;
use crate::mp::ct::{CtEqual, CtOption, CtSelect};

/// A concrete ECDSA instantiation: a specific curve paired with a specific
/// message hash function, plus the one curve operation ECDSA needs on top of
/// the generic group law.
///
/// This is implemented by dataless marker types — one per (curve, hash)
/// pair — and bundles:
///
/// * the [`Point`](Self::Point) curve group the scheme operates on (a
///   [`CurveGroup`]) and its [`Scalar`](Self::Scalar) field;
/// * [`x_mod_n`](Self::x_mod_n), the affine x-coordinate of a point reduced
///   modulo the group order;
/// * [`hash_to_scalar`](Self::hash_to_scalar), hashing a message into a scalar.
///
/// The provided marker types cover the SHA-2 family on every enabled curve
/// (e.g. [`P256R1_Sha256`] for SHA-256 with NIST P-256); see the module
/// documentation for plugging in a custom hash function.
pub trait EcdsaOperations {
    /// The curve group (point type) the scheme operates on.
    type Point: CurveGroup<Scalar = Self::Scalar>;

    /// The scalar field of the curve: the field a message hashes into and
    /// that the signature components live in. Always equal to the scalar
    /// field of [`Point`](Self::Point).
    type Scalar: Field;

    /// Hash a message into a scalar of the curve.
    ///
    /// Digests the message with this scheme's hash function and converts the
    /// digest to a scalar with the curve's `bits2int` transform.
    fn hash_to_scalar(message: &[u8]) -> Self::Scalar;

    /// The affine x-coordinate of `point` reduced modulo the group order, as a
    /// [`CtOption`] marked not-present for the identity (the point at
    /// infinity).
    ///
    /// The conversion out of the point representation must be constant time:
    /// during signing the point is secret-derived, and in particular whether
    /// it is the identity must not leak through a branch. The final
    /// reduction of the x-coordinate may be variable time, as that value is
    /// revealed as the `r` component of the signature anyway.
    fn x_mod_n(point: &Self::Point) -> CtOption<Self::Scalar>;
}

/// An ECDSA signature: the pair of scalars `(r, s)`, both non-zero.
///
/// Fixed size serialization (`to_bytes` / `from_bytes`) is provided for each
/// curve's scalar type.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Signature<S> {
    r: S,
    s: S,
}

impl<S: Field> Signature<S> {
    /// Create a signature from its two components.
    ///
    /// Returns `None` if either component is zero, as no valid signature has
    /// a zero component.
    pub fn from_scalars(r: S, s: S) -> Option<Self> {
        if r.is_zero() || s.is_zero() {
            return None;
        }
        Some(Signature { r, s })
    }

    /// The `r` component
    pub fn r(&self) -> &S {
        &self.r
    }

    /// The `s` component
    pub fn s(&self) -> &S {
        &self.s
    }
}

/// Compute the public key point associated with a secret key
pub fn public_key<G: CurveGroup>(secret: &G::Scalar) -> G {
    G::mul_base(secret)
}

/// Sign a message digest already converted to a scalar.
///
/// `nonce` must be unique and uniformly random for every signature; see the
/// module documentation.
///
/// The computation is constant time with respect to `secret`, `nonce` and
/// `hashed`: nothing branches on them, and the failure cases are folded into
/// the presence choice of the returned [`CtOption`] instead. The one
/// exception is the reduction of the x-coordinate of `nonce * G` into `r`,
/// which may be variable time on the value it produces, the public `r`
/// component of the signature (see [`EcdsaOperations::x_mod_n`]). The result
/// is marked not-present when `secret` or `nonce` is zero, or for the
/// negligibly rare nonces leading to a zero `r` or `s` component (retry with
/// a fresh nonce).
pub fn sign_hashed<O: EcdsaOperations>(
    secret: &O::Scalar,
    nonce: &O::Scalar,
    hashed: O::Scalar,
) -> CtOption<Signature<O::Scalar>> {
    let zero = &O::Scalar::ZERO;
    let nonce_nonzero = nonce.ct_ne(zero);

    // R = k*G ; r = R.x mod n (r_present is false only for a zero nonce)
    let (r_present, r) = O::x_mod_n(&O::Point::mul_base(nonce)).into_parts();

    // s = k^-1 * (z + r*d) mod n ; a zero nonce, whose result is discarded
    // through the validity choice, is substituted with 1 so that the
    // inversion stays defined
    let k = O::Scalar::ct_select(nonce_nonzero, nonce, &O::Scalar::ONE);
    let s = k.inverse() * (hashed + r.clone() * secret);

    let valid = secret.ct_ne(zero) & nonce_nonzero & r_present & r.ct_ne(zero) & s.ct_ne(zero);
    CtOption::from((valid, Signature { r, s }))
}

/// Sign a message: hash it into a scalar with `O`'s hash function, then sign
/// the resulting scalar.
///
/// `nonce` must be unique and uniformly random for every signature; see the
/// module documentation. Signing is constant time and reports failure through
/// the returned [`CtOption`]; see [`sign_hashed`].
pub fn sign<O: EcdsaOperations>(
    secret: &O::Scalar,
    nonce: &O::Scalar,
    message: &[u8],
) -> CtOption<Signature<O::Scalar>> {
    sign_hashed::<O>(secret, nonce, O::hash_to_scalar(message))
}

/// Verify a signature over a message digest already converted to a scalar.
///
/// `public` is assumed to be a valid public key: a point on the curve
/// (guaranteed by construction of the point type) that is not the point at
/// infinity.
pub fn verify_hashed<O: EcdsaOperations>(
    public: &O::Point,
    hashed: O::Scalar,
    signature: &Signature<O::Scalar>,
) -> bool {
    // r and s are non-zero canonical scalars by construction of Signature
    // R = (z/s)*G + (r/s)*Q ; accept only if R.x mod n == r
    let s_inv = signature.s.inverse();
    let u1 = hashed * &s_inv;
    let u2 = signature.r.clone() * s_inv;
    let rp = O::Point::mul_base(&u1) + public.mul_vartime(&u2);
    // verification only handles public data, so it can leave the
    // constant-time domain
    match O::x_mod_n(&rp).into_option() {
        None => false, // the identity has no x-coordinate
        Some(x) => x == signature.r,
    }
}

/// Verify a signature over a message: hash it into a scalar with `O`'s hash
/// function, then verify against the resulting scalar.
///
/// `public` is assumed to be a valid public key; see [`verify_hashed`].
pub fn verify<O: EcdsaOperations>(
    public: &O::Point,
    message: &[u8],
    signature: &Signature<O::Scalar>,
) -> bool {
    verify_hashed::<O>(public, O::hash_to_scalar(message), signature)
}

// the sec2 curve features that ecdsa can be instantiated on
macro_rules! cfg_any_curve {
    ($($item:item)*) => {
        $(
            #[cfg(any(
                feature = "p192k1",
                feature = "p192r1",
                feature = "p224k1",
                feature = "p224r1",
                feature = "p256k1",
                feature = "p256r1",
                feature = "p384r1",
                feature = "p521r1",
            ))]
            $item
        )*
    };
}

cfg_any_curve! {

/// Shift a big-endian byte buffer right by `bits` bits (0 to 7), in place.
fn shr_be(buf: &mut [u8], bits: usize) {
    debug_assert!(bits < 8);
    if bits == 0 {
        return;
    }
    let mut carry = 0u8;
    for b in buf.iter_mut() {
        let new_carry = *b << (8 - bits);
        *b = (*b >> bits) | carry;
        carry = new_carry;
    }
}

/// Generate the [`EcdsaOperations`] marker types pairing one curve with the
/// given cryptoxide hash functions.
macro_rules! ecdsa_hashing_impl {
    ($curve:ident, $($name:ident => $hasher:ident),+ $(,)?) => {
        $(
            #[allow(non_camel_case_types)]
            #[doc = concat!(
                "[`EcdsaOperations`] marker: ECDSA on [`", stringify!($curve),
                "`](crate::curve::sec2::", stringify!($curve),
                ") with ", stringify!($hasher), " message hashing"
            )]
            pub struct $name;

            impl EcdsaOperations for $name {
                type Point = crate::curve::sec2::$curve::Point;
                type Scalar = crate::curve::sec2::$curve::Scalar;

                fn hash_to_scalar(message: &[u8]) -> Self::Scalar {
                    $curve::digest_to_scalar(
                        &cryptoxide::hashing::sha2::$hasher::new().update(message).finalize(),
                    )
                }

                fn x_mod_n(point: &Self::Point) -> crate::mp::ct::CtOption<Self::Scalar> {
                    $curve::x_mod_n(point)
                }
            }
        )+
    };
}

/// Generate the per-curve items backing the curve's [`EcdsaOperations`]
/// markers: the `x_mod_n` point conversion, the bits2int `digest_to_scalar`
/// conversion, and the fixed size signature serialization.
macro_rules! ecdsa_impl {
    ($curve:path) => {
        use $curve::{FieldElement, Point, Scalar};

        /// Size in bytes of a serialized signature (`r || s`)
        pub const SIGNATURE_SIZE: usize = Scalar::SIZE_BYTES * 2;

        /// Reduce a scalar-sized big-endian buffer modulo the curve order.
        ///
        /// The value is at most `SIZE_BITS` bits, so it is either already
        /// canonical — covered by the cheap decode — or a single reduction
        /// away from it: only the rare values at or above the order (at most
        /// a `(2^SIZE_BITS - n) / 2^SIZE_BITS` fraction) pay for the wide
        /// byte-per-byte reduction. The variable timing is fine as every
        /// value passed here is public: a message digest, or an x-coordinate
        /// revealed as the `r` component of the signature.
        fn reduce_bytes_be(buf: [u8; Scalar::SIZE_BYTES]) -> Scalar {
            match Scalar::from_slice_be(&buf) {
                Some(s) => s,
                None => {
                    let mut wide = [0u8; Scalar::SIZE_BYTES * 2];
                    wide[Scalar::SIZE_BYTES..].copy_from_slice(&buf);
                    Scalar::init_from_wide_bytes_be(wide)
                }
            }
        }

        /// Convert a digest to a scalar with the SEC1 `bits2int` transform:
        /// keep the leftmost `min(qlen, 8*N)` bits of the digest interpreted
        /// as a big-endian integer, then reduce modulo the curve order.
        ///
        /// The branch is on constants, so each digest size monomorphizes to
        /// a single path; both end in [`reduce_bytes_be`], whose reduction
        /// is skipped entirely when the digest is strictly shorter than the
        /// order (its value is then always canonical).
        pub fn digest_to_scalar<const N: usize>(digest: &[u8; N]) -> Scalar {
            let mut buf = [0u8; Scalar::SIZE_BYTES];
            if N * 8 <= Scalar::SIZE_BITS {
                // the digest fits whole: left-pad to the scalar size
                buf[Scalar::SIZE_BYTES - N..].copy_from_slice(digest);
            } else {
                // keep the leftmost SIZE_BITS bits: the first SIZE_BYTES bytes
                // shifted right by the number of bits SIZE_BYTES overshoots
                buf.copy_from_slice(&digest[..Scalar::SIZE_BYTES]);
                super::shr_be(&mut buf, Scalar::SIZE_BYTES * 8 - Scalar::SIZE_BITS);
            }
            reduce_bytes_be(buf)
        }

        /// Reduce an x-coordinate modulo the curve order.
        ///
        /// The field and the scalar field can have different byte sizes, and
        /// the constant branch picks the one conversion the curve needs: a
        /// field no bigger than the scalar left-pads and reduces (when it is
        /// strictly smaller, e.g. p224k1 whose order is 225 bits, the value
        /// is always canonical and only the cheap decode runs); the wide
        /// reduction remains for a hypothetical field wider than the order,
        /// which no enabled curve has.
        fn field_to_scalar(fe: &FieldElement) -> Scalar {
            if FieldElement::SIZE_BYTES <= Scalar::SIZE_BYTES {
                let mut buf = [0u8; Scalar::SIZE_BYTES];
                buf[Scalar::SIZE_BYTES - FieldElement::SIZE_BYTES..]
                    .copy_from_slice(&fe.to_bytes_be());
                reduce_bytes_be(buf)
            } else {
                let mut wide = [0u8; Scalar::SIZE_BYTES * 2];
                let wlen = wide.len();
                wide[wlen - FieldElement::SIZE_BYTES..].copy_from_slice(&fe.to_bytes_be());
                Scalar::init_from_wide_bytes_be(wide)
            }
        }

        /// The affine x-coordinate of `point` reduced modulo the group order,
        /// as a [`CtOption`](crate::mp::ct::CtOption) marked not-present for
        /// the identity (the point at infinity). Backs
        /// [`EcdsaOperations::x_mod_n`](crate::protocol::ecdsa::EcdsaOperations::x_mod_n);
        /// see there for the constant-time contract.
        pub fn x_mod_n(point: &Point) -> crate::mp::ct::CtOption<Scalar> {
            let (present, x) = point.to_affine_x_ct().into_parts();
            crate::mp::ct::CtOption::from((present, field_to_scalar(&x)))
        }

        impl crate::protocol::ecdsa::Signature<Scalar> {
            /// Serialize to the fixed size `r || s` big-endian representation
            pub fn to_bytes(&self) -> [u8; SIGNATURE_SIZE] {
                let mut out = [0u8; SIGNATURE_SIZE];
                out[..Scalar::SIZE_BYTES].copy_from_slice(&self.r().to_bytes_be());
                out[Scalar::SIZE_BYTES..].copy_from_slice(&self.s().to_bytes_be());
                out
            }

            /// Deserialize from the fixed size `r || s` big-endian
            /// representation, rejecting non-canonical (`>= order`) or zero
            /// components
            pub fn from_bytes(bytes: &[u8; SIGNATURE_SIZE]) -> Option<Self> {
                let r = Scalar::from_slice_be(&bytes[..Scalar::SIZE_BYTES])?;
                let s = Scalar::from_slice_be(&bytes[Scalar::SIZE_BYTES..])?;
                Self::from_scalars(r, s)
            }
        }
    };
}

}

#[cfg(feature = "p192k1")]
/// ECDSA support for the [`crate::curve::sec2::p192k1`] curve
pub mod p192k1 {
    ecdsa_impl!(crate::curve::sec2::p192k1);
}
#[cfg(feature = "p192k1")]
ecdsa_hashing_impl!(
    p192k1,
    P192K1_Sha224 => Sha224,
    P192K1_Sha256 => Sha256,
    P192K1_Sha384 => Sha384,
    P192K1_Sha512 => Sha512,
);

#[cfg(feature = "p192r1")]
/// ECDSA support for the [`crate::curve::sec2::p192r1`] curve (NIST P-192)
pub mod p192r1 {
    ecdsa_impl!(crate::curve::sec2::p192r1);
}
#[cfg(feature = "p192r1")]
ecdsa_hashing_impl!(
    p192r1,
    P192R1_Sha224 => Sha224,
    P192R1_Sha256 => Sha256,
    P192R1_Sha384 => Sha384,
    P192R1_Sha512 => Sha512,
);

#[cfg(feature = "p224k1")]
/// ECDSA support for the [`crate::curve::sec2::p224k1`] curve
pub mod p224k1 {
    ecdsa_impl!(crate::curve::sec2::p224k1);
}
#[cfg(feature = "p224k1")]
ecdsa_hashing_impl!(
    p224k1,
    P224K1_Sha224 => Sha224,
    P224K1_Sha256 => Sha256,
    P224K1_Sha384 => Sha384,
    P224K1_Sha512 => Sha512,
);

#[cfg(feature = "p224r1")]
/// ECDSA support for the [`crate::curve::sec2::p224r1`] curve (NIST P-224)
pub mod p224r1 {
    ecdsa_impl!(crate::curve::sec2::p224r1);
}
#[cfg(feature = "p224r1")]
ecdsa_hashing_impl!(
    p224r1,
    P224R1_Sha224 => Sha224,
    P224R1_Sha256 => Sha256,
    P224R1_Sha384 => Sha384,
    P224R1_Sha512 => Sha512,
);

#[cfg(feature = "p256k1")]
/// ECDSA support for the [`crate::curve::sec2::p256k1`] curve (bitcoin's secp256k1)
pub mod p256k1 {
    ecdsa_impl!(crate::curve::sec2::p256k1);
}
#[cfg(feature = "p256k1")]
ecdsa_hashing_impl!(
    p256k1,
    P256K1_Sha224 => Sha224,
    P256K1_Sha256 => Sha256,
    P256K1_Sha384 => Sha384,
    P256K1_Sha512 => Sha512,
);

#[cfg(feature = "p256r1")]
/// ECDSA support for the [`crate::curve::sec2::p256r1`] curve (NIST P-256)
pub mod p256r1 {
    ecdsa_impl!(crate::curve::sec2::p256r1);
}
#[cfg(feature = "p256r1")]
ecdsa_hashing_impl!(
    p256r1,
    P256R1_Sha224 => Sha224,
    P256R1_Sha256 => Sha256,
    P256R1_Sha384 => Sha384,
    P256R1_Sha512 => Sha512,
);

#[cfg(feature = "p384r1")]
/// ECDSA support for the [`crate::curve::sec2::p384r1`] curve (NIST P-384)
pub mod p384r1 {
    ecdsa_impl!(crate::curve::sec2::p384r1);
}
#[cfg(feature = "p384r1")]
ecdsa_hashing_impl!(
    p384r1,
    P384R1_Sha224 => Sha224,
    P384R1_Sha256 => Sha256,
    P384R1_Sha384 => Sha384,
    P384R1_Sha512 => Sha512,
);

#[cfg(feature = "p521r1")]
/// ECDSA support for the [`crate::curve::sec2::p521r1`] curve (NIST P-521)
pub mod p521r1 {
    ecdsa_impl!(crate::curve::sec2::p521r1);
}
#[cfg(feature = "p521r1")]
ecdsa_hashing_impl!(
    p521r1,
    P521R1_Sha224 => Sha224,
    P521R1_Sha256 => Sha256,
    P521R1_Sha384 => Sha384,
    P521R1_Sha512 => Sha512,
);

#[cfg(test)]
mod tests {
    fn hex(s: &str, size: usize) -> Vec<u8> {
        let s = if s.len() % 2 == 1 {
            format!("0{}", s)
        } else {
            s.to_string()
        };
        let v: Vec<u8> = (0..s.len() / 2)
            .map(|i| u8::from_str_radix(&s[i * 2..i * 2 + 2], 16).unwrap())
            .collect();
        assert!(v.len() <= size);
        let mut out = vec![0u8; size - v.len()];
        out.extend(v);
        out
    }

    #[derive(Clone, Copy)]
    enum Alg {
        Sha224,
        Sha256,
        Sha384,
        Sha512,
    }

    struct Kat {
        alg: Alg,
        message: &'static str,
        k: &'static str,
        r: &'static str,
        s: &'static str,
    }

    #[test]
    fn shr_be() {
        let mut buf = [0x81, 0xff, 0x00, 0x01];
        super::shr_be(&mut buf, 0);
        assert_eq!(buf, [0x81, 0xff, 0x00, 0x01]);
        super::shr_be(&mut buf, 1);
        assert_eq!(buf, [0x40, 0xff, 0x80, 0x00]);
        let mut buf = [0x81, 0xff, 0x00, 0x01];
        super::shr_be(&mut buf, 7);
        assert_eq!(buf, [0x01, 0x03, 0xfe, 0x00]);
    }

    /// Sign/verify roundtrip, tamper rejection and serialization tests for
    /// one curve, generated inside that curve's test module (which aliases
    /// its hashing markers to `Sha256` / `Sha512`).
    macro_rules! ecdsa_test_roundtrip {
        () => {
            #[test]
            fn roundtrip_and_tamper() {
                let secret = Scalar::init_from_wide_bytes_be([0x42; Scalar::SIZE_BYTES * 2]);
                let nonce = Scalar::init_from_wide_bytes_be([0xac; Scalar::SIZE_BYTES * 2]);
                let public: Point = public_key(&secret);
                let msg = b"attack at dawn";

                // exercise both the pad (sha256 on larger curves) and
                // truncate (sha512) paths of digest_to_scalar
                let sig = sign::<Sha256>(&secret, &nonce, msg).into_option().unwrap();
                assert!(verify::<Sha256>(&public, msg, &sig));
                let sig512 = sign::<Sha512>(&secret, &nonce, msg).into_option().unwrap();
                assert!(verify::<Sha512>(&public, msg, &sig512));

                // wrong message / wrong hash / wrong key rejected
                assert!(!verify::<Sha256>(&public, b"attack at dusk", &sig));
                assert!(!verify::<Sha256>(&public, msg, &sig512));
                let other: Point = public_key(&Scalar::from_u64(1234));
                assert!(!verify::<Sha256>(&other, msg, &sig));

                // the hashing marker agrees with hashing manually, converting
                // with digest_to_scalar and using the scalar-level functions
                let digest = cryptoxide::hashing::sha2::Sha256::new()
                    .update(msg)
                    .finalize();
                let z = digest_to_scalar(&digest);
                assert_eq!(
                    sign_hashed::<Sha256>(&secret, &nonce, z.clone())
                        .into_option()
                        .unwrap(),
                    sig
                );
                assert!(verify_hashed::<Sha256>(&public, z, &sig));

                // serialization roundtrip
                let restored = Signature::<Scalar>::from_bytes(&sig.to_bytes()).unwrap();
                assert_eq!(restored, sig);
                assert!(verify::<Sha256>(&public, msg, &restored));

                // zero components rejected
                assert!(Signature::<Scalar>::from_bytes(&[0u8; SIGNATURE_SIZE]).is_none());
                let mut half_zero = sig.to_bytes();
                for b in half_zero[Scalar::SIZE_BYTES..].iter_mut() {
                    *b = 0
                }
                assert!(Signature::<Scalar>::from_bytes(&half_zero).is_none());

                // zero secret or nonce rejected
                assert!(sign::<Sha256>(&Scalar::zero(), &nonce, msg)
                    .into_option()
                    .is_none());
                assert!(sign::<Sha256>(&secret, &Scalar::zero(), msg)
                    .into_option()
                    .is_none());
            }

            /// Check the fast canonical-decode path and its rare wide
            /// fallback in digest_to_scalar against an always-wide Horner
            /// reduction reference
            #[test]
            fn digest_reduction_reference() {
                let c256 = Scalar::from_u64(256);
                let horner = |bytes: &[u8]| {
                    let mut acc = Scalar::zero();
                    for b in bytes {
                        acc = &acc * &c256 + Scalar::from_u64(*b as u64);
                    }
                    acc
                };
                // bits2int truncation, then the reference reduction
                let reference = |digest: &[u8]| {
                    if digest.len() * 8 <= Scalar::SIZE_BITS {
                        horner(digest)
                    } else {
                        let mut buf = [0u8; Scalar::SIZE_BYTES];
                        buf.copy_from_slice(&digest[..Scalar::SIZE_BYTES]);
                        crate::protocol::ecdsa::shr_be(
                            &mut buf,
                            Scalar::SIZE_BYTES * 8 - Scalar::SIZE_BITS,
                        );
                        horner(&buf)
                    }
                };

                // short digest: strictly below the order, canonical decode only
                assert_eq!(digest_to_scalar(&[0xfeu8; 8]), reference(&[0xfeu8; 8]));
                // scalar-sized all-ones digest: maximal value, hitting the
                // non-canonical fallback on the byte-aligned curves and the
                // truncation path on the others (225/521-bit orders)
                let ones = [0xffu8; Scalar::SIZE_BYTES];
                assert_eq!(digest_to_scalar(&ones), reference(&ones));
                // double-sized all-ones digest: truncation on every curve
                let wide_ones = [0xffu8; Scalar::SIZE_BYTES * 2];
                assert_eq!(digest_to_scalar(&wide_ones), reference(&wide_ones));
            }
        };
    }

    /// KAT + roundtrip tests for one curve.
    ///
    /// The KAT vectors come from RFC 6979: even though the nonce is a
    /// parameter here (this module does not do RFC 6979 nonce derivation),
    /// the RFC publishes the derived `k` for each vector, so feeding that `k`
    /// back must reproduce the published `(r, s)` exactly.
    macro_rules! ecdsa_test {
        ($module:ident, $sha224:ident, $sha256:ident, $sha384:ident, $sha512:ident,
         $secret:expr, $ux:expr, $uy:expr, $kats:expr) => {
            mod $module {
                use super::super::$module::{digest_to_scalar, SIGNATURE_SIZE};
                use super::{hex, Alg, Kat};
                use crate::curve::sec2::$module::{FieldElement, Point, Scalar};
                use crate::protocol::ecdsa::{
                    public_key, sign, sign_hashed, verify, verify_hashed, $sha224 as Sha224,
                    $sha256 as Sha256, $sha384 as Sha384, $sha512 as Sha512, Signature,
                };

                fn scalar(s: &str) -> Scalar {
                    Scalar::from_slice_be(&hex(s, Scalar::SIZE_BYTES)).unwrap()
                }

                fn fe(s: &str) -> FieldElement {
                    FieldElement::from_slice_be(&hex(s, FieldElement::SIZE_BYTES)).unwrap()
                }

                ecdsa_test_roundtrip!();

                #[test]
                fn rfc6979_kats() {
                    let secret = scalar($secret);
                    let public: Point = public_key(&secret);

                    // the RFC also publishes the public key point
                    let pa = public.to_affine().unwrap();
                    let (x, y) = pa.to_coordinate();
                    assert_eq!(x, &fe($ux));
                    assert_eq!(y, &fe($uy));

                    for kat in $kats {
                        let k = scalar(kat.k);
                        let msg = kat.message.as_bytes();
                        let sig = match kat.alg {
                            Alg::Sha224 => sign::<Sha224>(&secret, &k, msg),
                            Alg::Sha256 => sign::<Sha256>(&secret, &k, msg),
                            Alg::Sha384 => sign::<Sha384>(&secret, &k, msg),
                            Alg::Sha512 => sign::<Sha512>(&secret, &k, msg),
                        }
                        .into_option()
                        .unwrap();
                        assert_eq!(sig.r(), &scalar(kat.r), "r mismatch {}", kat.message);
                        assert_eq!(sig.s(), &scalar(kat.s), "s mismatch {}", kat.message);

                        let ok = match kat.alg {
                            Alg::Sha224 => verify::<Sha224>(&public, msg, &sig),
                            Alg::Sha256 => verify::<Sha256>(&public, msg, &sig),
                            Alg::Sha384 => verify::<Sha384>(&public, msg, &sig),
                            Alg::Sha512 => verify::<Sha512>(&public, msg, &sig),
                        };
                        assert!(ok, "valid signature rejected {}", kat.message);
                    }
                }
            }
        };
    }

    // RFC 6979 A.2.3, NIST P-192
    #[cfg(feature = "p192r1")]
    ecdsa_test!(
        p192r1,
        P192R1_Sha224,
        P192R1_Sha256,
        P192R1_Sha384,
        P192R1_Sha512,
        "6FAB034934E4C0FC9AE67F5B5659A9D7D1FEFD187EE09FD4",
        "AC2C77F529F91689FEA0EA5EFEC7F210D8EEA0B9E047ED56",
        "3BC723E57670BD4887EBC732C523063D0A7C957BC97C1C43",
        &[
            Kat {
                alg: Alg::Sha256,
                message: "sample",
                k: "32B1B6D7D42A05CB449065727A84804FB1A3E34D8F261496",
                r: "4B0B8CE98A92866A2820E20AA6B75B56382E0F9BFD5ECB55",
                s: "CCDB006926EA9565CBADC840829D8C384E06DE1F1E381B85",
            },
            Kat {
                alg: Alg::Sha512,
                message: "sample",
                k: "A2AC7AB055E4F20692D49209544C203A7D1F2C0BFBC75DB1",
                r: "4D60C5AB1996BD848343B31C00850205E2EA6922DAC2E4B8",
                s: "3F6E837448F027A1BF4B34E796E32A811CBB4050908D8F67",
            },
            Kat {
                alg: Alg::Sha256,
                message: "test",
                k: "5C4CE89CF56D9E7C77C8585339B006B97B5F0680B4306C6C",
                r: "3A718BD8B4926C3B52EE6BBE67EF79B18CB6EB62B1AD97AE",
                s: "5662E6848A4A19B1F1AE2F72ACD4B8BBE50F1EAC65D9124F",
            },
        ]
    );

    // RFC 6979 A.2.4, NIST P-224
    #[cfg(feature = "p224r1")]
    ecdsa_test!(
        p224r1,
        P224R1_Sha224,
        P224R1_Sha256,
        P224R1_Sha384,
        P224R1_Sha512,
        "F220266E1105BFE3083E03EC7A3A654651F45E37167E88600BF257C1",
        "00CF08DA5AD719E42707FA431292DEA11244D64FC51610D94B130D6C",
        "EEAB6F3DEBE455E3DBF85416F7030CBD94F34F2D6F232C69F3C1385A",
        &[
            Kat {
                alg: Alg::Sha224,
                message: "sample",
                k: "C1D1F2F10881088301880506805FEB4825FE09ACB6816C36991AA06D",
                r: "1CDFE6662DDE1E4A1EC4CDEDF6A1F5A2FB7FBD9145C12113E6ABFD3E",
                s: "A6694FD7718A21053F225D3F46197CA699D45006C06F871808F43EBC",
            },
            Kat {
                alg: Alg::Sha256,
                message: "sample",
                k: "AD3029E0278F80643DE33917CE6908C70A8FF50A411F06E41DEDFCDC",
                r: "61AA3DA010E8E8406C656BC477A7A7189895E7E840CDFE8FF42307BA",
                s: "BC814050DAB5D23770879494F9E0A680DC1AF7161991BDE692B10101",
            },
            Kat {
                alg: Alg::Sha512,
                message: "test",
                k: "E39C2AA4EA6BE2306C72126D40ED77BF9739BB4D6EF2BBB1DCB6169D",
                r: "049F050477C5ADD858CAC56208394B5A55BAEBBE887FDF765047C17C",
                s: "077EB13E7005929CEFA3CD0403C7CDCC077ADF4E44F3C41B2F60ECFF",
            },
        ]
    );

    // RFC 6979 A.2.5, NIST P-256
    #[cfg(feature = "p256r1")]
    ecdsa_test!(
        p256r1,
        P256R1_Sha224,
        P256R1_Sha256,
        P256R1_Sha384,
        P256R1_Sha512,
        "C9AFA9D845BA75166B5C215767B1D6934E50C3DB36E89B127B8A622B120F6721",
        "60FED4BA255A9D31C961EB74C6356D68C049B8923B61FA6CE669622E60F29FB6",
        "7903FE1008B8BC99A41AE9E95628BC64F2F1B20C2D7E9F5177A3C294D4462299",
        &[
            Kat {
                alg: Alg::Sha256,
                message: "sample",
                k: "A6E3C57DD01ABE90086538398355DD4C3B17AA873382B0F24D6129493D8AAD60",
                r: "EFD48B2AACB6A8FD1140DD9CD45E81D69D2C877B56AAF991C34D0EA84EAF3716",
                s: "F7CB1C942D657C41D436C7A1B6E29F65F3E900DBB9AFF4064DC4AB2F843ACDA8",
            },
            Kat {
                alg: Alg::Sha512,
                message: "sample",
                k: "5FA81C63109BADB88C1F367B47DA606DA28CAD69AA22C4FE6AD7DF73A7173AA5",
                r: "8496A60B5E9B47C825488827E0495B0E3FA109EC4568FD3F8D1097678EB97F00",
                s: "2362AB1ADBE2B8ADF9CB9EDAB740EA6049C028114F2460F96554F61FAE3302FE",
            },
            Kat {
                alg: Alg::Sha256,
                message: "test",
                k: "D16B6AE827F17175E040871A1C7EC3500192C4C92677336EC2537ACAEE0008E0",
                r: "F1ABB023518351CD71D881567B1EA663ED3EFCF6C5132B354F28D3B0B7D38367",
                s: "019F4113742A2B14BD25926B49C649155F267E60D3814B4C0CC84250E46F0083",
            },
        ]
    );

    // RFC 6979 A.2.6, NIST P-384
    #[cfg(feature = "p384r1")]
    ecdsa_test!(
        p384r1,
        P384R1_Sha224,
        P384R1_Sha256,
        P384R1_Sha384,
        P384R1_Sha512,
        "6B9D3DAD2E1B8C1C05B19875B6659F4DE23C3B667BF297BA9AA47740787137D896D5724E4C70A825F872C9EA60D2EDF5",
        "EC3A4E415B4E19A4568618029F427FA5DA9A8BC4AE92E02E06AAE5286B300C64DEF8F0EA9055866064A254515480BC13",
        "8015D9B72D7D57244EA8EF9AC0C621896708A59367F9DFB9F54CA84B3F1C9DB1288B231C3AE0D4FE7344FD2533264720",
        &[
            Kat {
                alg: Alg::Sha256,
                message: "sample",
                k: "180AE9F9AEC5438A44BC159A1FCB277C7BE54FA20E7CF404B490650A8ACC414E375572342863C899F9F2EDF9747A9B60",
                r: "21B13D1E013C7FA1392D03C5F99AF8B30C570C6F98D4EA8E354B63A21D3DAA33BDE1E888E63355D92FA2B3C36D8FB2CD",
                s: "F3AA443FB107745BF4BD77CB3891674632068A10CA67E3D45DB2266FA7D1FEEBEFDC63ECCD1AC42EC0CB8668A4FA0AB0",
            },
            Kat {
                alg: Alg::Sha384,
                message: "sample",
                k: "94ED910D1A099DAD3254E9242AE85ABDE4BA15168EAF0CA87A555FD56D10FBCA2907E3E83BA95368623B8C4686915CF9",
                r: "94EDBB92A5ECB8AAD4736E56C691916B3F88140666CE9FA73D64C4EA95AD133C81A648152E44ACF96E36DD1E80FABE46",
                s: "99EF4AEB15F178CEA1FE40DB2603138F130E740A19624526203B6351D0A3A94FA329C145786E679E7B82C71A38628AC8",
            },
            Kat {
                alg: Alg::Sha512,
                message: "test",
                k: "3780C4F67CB15518B6ACAE34C9F83568D2E12E47DEAB6C50A4E4EE5319D1E8CE0E2CC8A136036DC4B9C00E6888F66B6C",
                r: "A0D5D090C9980FAF3C2CE57B7AE951D31977DD11C775D314AF55F76C676447D06FB6495CD21B4B6E340FC236584FB277",
                s: "976984E59B4C77B0E8E4460DCA3D9F20E07B9BB1F63BEEFAF576F6B2E8B224634A2092CD3792E0159AD9CEE37659C736",
            },
        ]
    );

    // RFC 6979 A.2.7, NIST P-521
    #[cfg(feature = "p521r1")]
    ecdsa_test!(
        p521r1,
        P521R1_Sha224,
        P521R1_Sha256,
        P521R1_Sha384,
        P521R1_Sha512,
        "0FAD06DAA62BA3B25D2FB40133DA757205DE67F5BB0018FEE8C86E1B68C7E75CAA896EB32F1F47C70855836A6D16FCC1466F6D8FBEC67DB89EC0C08B0E996B83538",
        "1894550D0785932E00EAA23B694F213F8C3121F86DC97A04E5A7167DB4E5BCD371123D46E45DB6B5D5370A7F20FB633155D38FFA16D2BD761DCAC474B9A2F5023A4",
        "0493101C962CD4D2FDDF782285E64584139C2F91B47F87FF82354D6630F746A28A0DB25741B5B34A828008B22ACC23F924FAAFBD4D33F81EA66956DFEAA2BFDFCF5",
        &[
            Kat {
                alg: Alg::Sha256,
                message: "sample",
                k: "0EDF38AFCAAECAB4383358B34D67C9F2216C8382AAEA44A3DAD5FDC9C32575761793FEF24EB0FC276DFC4F6E3EC476752F043CF01415387470BCBD8678ED2C7E1A0",
                r: "1511BB4D675114FE266FC4372B87682BAECC01D3CC62CF2303C92B3526012659D16876E25C7C1E57648F23B73564D67F61C6F14D527D54972810421E7D87589E1A7",
                s: "04A171143A83163D6DF460AAF61522695F207A58B95C0644D87E52AA1A347916E4F7A72930B1BC06DBE22CE3F58264AFD23704CBB63B29B931F7DE6C9D949A7ECFC",
            },
            Kat {
                alg: Alg::Sha512,
                message: "sample",
                k: "1DAE2EA071F8110DC26882D4D5EAE0621A3256FC8847FB9022E2B7D28E6F10198B1574FDD03A9053C08A1854A168AA5A57470EC97DD5CE090124EF52A2F7ECBFFD3",
                r: "0C328FAFCBD79DD77850370C46325D987CB525569FB63C5D3BC53950E6D4C5F174E25A1EE9017B5D450606ADD152B534931D7D4E8455CC91F9B15BF05EC36E377FA",
                s: "0617CCE7CF5064806C467F678D3B4080D6F1CC50AF26CA209417308281B68AF282623EAA63E5B5C0723D8B8C37FF0777B1A20F8CCB1DCCC43997F1EE0E44DA4A67A",
            },
        ]
    );

    // no standardized test vectors for the k1 curves; the shared macro code
    // is covered by the NIST KATs above, and these still get the full
    // roundtrip/tamper coverage (p224k1 notably exercises the non-byte-aligned
    // 225-bit order in both digest_to_scalar and field_to_scalar)
    macro_rules! ecdsa_test_roundtrip_only {
        ($module:ident, $sha256:ident, $sha512:ident) => {
            mod $module {
                use super::super::$module::{digest_to_scalar, SIGNATURE_SIZE};
                use crate::curve::sec2::$module::{Point, Scalar};
                use crate::protocol::ecdsa::{
                    public_key, sign, sign_hashed, verify, verify_hashed, $sha256 as Sha256,
                    $sha512 as Sha512, Signature,
                };

                ecdsa_test_roundtrip!();
            }
        };
    }

    #[cfg(feature = "p192k1")]
    ecdsa_test_roundtrip_only!(p192k1, P192K1_Sha256, P192K1_Sha512);
    #[cfg(feature = "p224k1")]
    ecdsa_test_roundtrip_only!(p224k1, P224K1_Sha256, P224K1_Sha512);
    #[cfg(feature = "p256k1")]
    ecdsa_test_roundtrip_only!(p256k1, P256K1_Sha256, P256K1_Sha512);

    // a caller-supplied instantiation, exercising what a custom
    // `EcdsaOperations` implementation is built from: the curve's point type
    // and `x_mod_n`, plus an arbitrary digest run through the curve's
    // `digest_to_scalar` bits2int conversion
    #[cfg(feature = "p256r1")]
    mod custom_hashing {
        use crate::curve::sec2::p256r1::{Point, Scalar};
        use crate::mp::ct::CtOption;
        use crate::protocol::ecdsa::{
            p256r1::digest_to_scalar, public_key, sign, verify, EcdsaOperations, P256R1_Sha256,
        };

        struct DoubleSha256;

        impl EcdsaOperations for DoubleSha256 {
            type Point = Point;
            type Scalar = Scalar;

            fn hash_to_scalar(message: &[u8]) -> Scalar {
                use cryptoxide::hashing::sha2::Sha256;
                let digest = Sha256::new().update(message).finalize();
                digest_to_scalar(&Sha256::new().update(&digest).finalize())
            }

            fn x_mod_n(point: &Point) -> CtOption<Scalar> {
                // reuse the curve's own reduction
                crate::protocol::ecdsa::p256r1::x_mod_n(point)
            }
        }

        #[test]
        fn roundtrip() {
            let secret = Scalar::init_from_wide_bytes_be([0x42; Scalar::SIZE_BYTES * 2]);
            let nonce = Scalar::init_from_wide_bytes_be([0xac; Scalar::SIZE_BYTES * 2]);
            let public: Point = public_key(&secret);
            let msg = b"attack at dawn";

            let sig = sign::<DoubleSha256>(&secret, &nonce, msg)
                .into_option()
                .unwrap();
            assert!(verify::<DoubleSha256>(&public, msg, &sig));
            assert!(!verify::<DoubleSha256>(&public, b"attack at dusk", &sig));
            // a different hashing choice does not cross-verify
            assert!(!verify::<P256R1_Sha256>(&public, msg, &sig));
        }
    }
}
