//! Ed25519 signature scheme (RFC 8032, PureEdDSA).
//!
//! Keys, signatures and encoded points use the RFC 8032 little-endian wire
//! format. The internal SHA-512 is provided by the `cryptoxide` crate.
//!
//! # Alternative hash function
//!
//! Ed25519 is defined by RFC 8032 on top of SHA-512, and that is what the plain
//! API of this module uses. The scheme itself only needs a 512-bit hash
//! function though, so every operation that need a hash function is also available
//! parametrized through the `_with` suffixed methods (e.g. [`Keypair::sign_with`]),
//! where the hash function is any `Fn(&[&[u8]]) -> [u8; HASH_LENGTH]`
//! hashing the concatenation of all the given slices.
//!
//! Two such hash functions are provided: [`sha512`] (the standard one) and,
//! with the `ed25519-blake2` feature, BLAKE2b-512.
#![cfg_attr(
    feature = "ed25519-blake2",
    doc = "For convenience, the [`blake2b`] module mirrors the API of this module pre-applied to [`blake2b512`]."
)]
//!
//! Note that keys and signatures are bound to the hash function they were
//! created with: they are *not* interchangeable with the standard SHA-512 ones.
//!
//! ```
//! use eccoxide::protocol::ed25519::{sha512, Keypair};
//!
//! // the standard scheme, spelled out through the parametrized interface
//! let keypair = Keypair::from_seed_with(sha512, [42u8; 32]);
//! let signature = keypair.sign_with(sha512, b"message");
//! assert!(keypair.public().verify(b"message", &signature));
//! ```

use crate::curve::curve25519::{FieldElement, Point, PrecomputedPoint, Scalar};
use crate::curve::field::Sign;
#[cfg(feature = "ed25519-blake2")]
use cryptoxide::hashing::blake2b::Blake2b;
use cryptoxide::hashing::sha2::Sha512;

/// Digest size (64 bytes) a hash function needs to output to be usable here
pub const HASH_LENGTH: usize = 64;

/// SHA-512 of the concatenation of all the given slices.
///
/// This is the hash function of the standard Ed25519 scheme, and thus the one
/// the plain API uses: handing it to the `_with` methods (e.g.
/// [`Keypair::sign_with`]) yields the same results as their plain counterpart.
pub fn sha512(slices: &[&[u8]]) -> [u8; HASH_LENGTH] {
    let mut h = Sha512::new();
    for p in slices {
        h.update_mut(p);
    }
    h.finalize()
}

/// BLAKE2b-512 of the concatenation of all the given slices.
///
/// A non standard choice of hash function for Ed25519, to be used with the
/// `_with` methods or through the [`blake2b`] module.
#[cfg(feature = "ed25519-blake2")]
pub fn blake2b512(slices: &[&[u8]]) -> [u8; HASH_LENGTH] {
    let mut h = Blake2b::<512>::new();
    for p in slices {
        h.update_mut(p);
    }
    h.finalize()
}

/// Reduce a 64-byte hash output, interpreted little-endian, modulo the group
/// order `l`.
fn reduce_wide_le(digest: &[u8; HASH_LENGTH]) -> Scalar {
    Scalar::init_from_wide_bytes_le(*digest)
}

/// Encode a point to its 32-byte RFC 8032 representation: little-endian y with
/// the low bit of x stored in the most significant bit.
fn encode_point(p: &Point) -> [u8; 32] {
    let (x, y) = p.to_affine();
    let mut out = y.to_bytes_le(); // little-endian, the native field order
    if let Sign::Negative = x.sign() {
        // Sign::Negative means the low bit of x is set
        out[31] |= 0x80;
    }
    out
}

/// Decode a 32-byte RFC 8032 point, rejecting non-canonical or off-curve inputs.
fn decode_point(bytes: &[u8; 32]) -> Option<Point> {
    let x_sign_bit = bytes[31] >> 7;
    let mut le = *bytes;
    le[31] &= 0x7f;
    // y must be a canonical field element (< p), little-endian on the wire
    let y = FieldElement::from_bytes_le(&le)?;
    // reject the (x = 0, sign = 1) non-canonical encoding. On the curve x is
    // zero exactly for y = ±1, since x = 0 forces y^2 = 1, so this is a test on
    // y — no need to recover x, which would cost a field inversion.
    if x_sign_bit == 1 {
        let one = FieldElement::one();
        if y == one || y == -&one {
            return None;
        }
    }
    let want = if x_sign_bit == 1 {
        Sign::Negative
    } else {
        Sign::Positive
    };
    Point::decompress(&y, want)
}

/// Expand a 32-byte seed into the secret scalar (mod l) and the nonce prefix,
/// with the given hash function.
fn expand_secret<H>(hash: &H, seed: &[u8; 32]) -> (Scalar, [u8; 32])
where
    H: Fn(&[&[u8]]) -> [u8; HASH_LENGTH],
{
    let h = hash(&[&seed[..]]);

    let mut a_le = [0u8; 32];
    a_le.copy_from_slice(&h[..32]);
    a_le[0] &= 248;
    a_le[31] &= 127;
    a_le[31] |= 64;

    // reduce the clamped little-endian scalar modulo l: place it in the low
    // half of the wide little-endian buffer (the high half stays zero)
    let mut wide = [0u8; 64];
    wide[..32].copy_from_slice(&a_le);
    let a = Scalar::init_from_wide_bytes_le(wide);

    let mut prefix = [0u8; 32];
    prefix.copy_from_slice(&h[32..]);
    (a, prefix)
}

fn public_from_seed<H>(hash: &H, seed: &[u8; 32]) -> [u8; 32]
where
    H: Fn(&[&[u8]]) -> [u8; HASH_LENGTH],
{
    let (a, _) = expand_secret(hash, seed);
    encode_point(&Point::mul_base(&a))
}

/// Core signing, given the expanded secret scalar `a`, the nonce `prefix`, and
/// the *already-encoded* public key `public` (A).
///
/// The signature only needs A inside the `k = H(R || A || M)` hash (A is not
/// part of the 64-byte `R || S` output), so when the public key is already known
/// — as it is for a [`Keypair`] — this avoids the extra fixed-base scalar
/// multiplication that recomputing A would cost.
fn sign_with_public<H>(
    hash: &H,
    a: &Scalar,
    prefix: &[u8; 32],
    public: &[u8; 32],
    message: &[u8],
) -> [u8; 64]
where
    H: Fn(&[&[u8]]) -> [u8; HASH_LENGTH],
{
    // r = H(prefix || M) mod l ; R = [r]B
    let r = reduce_wide_le(&hash(&[&prefix[..], message]));
    let r_encoded = encode_point(&Point::mul_base(&r));

    // k = H(R || A || M) mod l
    let k = reduce_wide_le(&hash(&[&r_encoded[..], &public[..], message]));

    // S = (r + k·a) mod l
    let s = r + &(&k * a);
    let s_le = s.to_bytes_le(); // little-endian, the native scalar order

    let mut sig = [0u8; 64];
    sig[..32].copy_from_slice(&r_encoded);
    sig[32..].copy_from_slice(&s_le);
    sig
}

fn sign<H>(hash: &H, seed: &[u8; 32], message: &[u8]) -> [u8; 64]
where
    H: Fn(&[&[u8]]) -> [u8; HASH_LENGTH],
{
    // A bare seed has no cached public key, so A must be derived here.
    let (a, prefix) = expand_secret(hash, seed);
    let public = encode_point(&Point::mul_base(&a));
    sign_with_public(hash, &a, &prefix, &public, message)
}

/// The verification equation common part for verify to cope with precomputed / non precomputed difference
///
/// `lhs` callback is handed `(S, k)` and must return `[S]B + [k](-A)`
fn verify_equation<H>(
    hash: &H,
    public: &[u8; 32],
    message: &[u8],
    sig: &[u8; 64],
    lhs: impl FnOnce(&Scalar, &Scalar) -> Point,
) -> bool
where
    H: Fn(&[&[u8]]) -> [u8; HASH_LENGTH],
{
    let mut r_encoded = [0u8; 32];
    r_encoded.copy_from_slice(&sig[..32]);
    let r_point = match decode_point(&r_encoded) {
        Some(p) => p,
        None => return false,
    };

    // S must be canonical (< l), little-endian on the wire
    let mut s_le = [0u8; 32];
    s_le.copy_from_slice(&sig[32..]);
    let s = match Scalar::from_bytes_le(&s_le) {
        Some(s) => s,
        None => return false,
    };

    // k = H(R || A || M) mod l
    let k = reduce_wide_le(&hash(&[&r_encoded[..], &public[..], message]));

    // accept iff [S]B == R + [k]A, rewritten as [S]B + [k](-A) == R so that
    // the two multiplications interleave into a variable op since all values are
    // public
    lhs(&s, &k) == r_point
}

fn verify<H>(hash: &H, public: &[u8; 32], message: &[u8], sig: &[u8; 64]) -> bool
where
    H: Fn(&[&[u8]]) -> [u8; HASH_LENGTH],
{
    let neg_a = match decode_point(public) {
        Some(p) => -&p,
        None => return false,
    };
    verify_equation(hash, public, message, sig, |s, k| {
        Point::double_scalar_mul_base_vartime(s, k, &neg_a)
    })
}

/// An Ed25519 secret key: the 32-byte seed from which everything is derived.
#[derive(Clone)]
pub struct SecretKey([u8; 32]);

/// An Ed25519 public key (32-byte compressed point).
#[derive(Clone, PartialEq, Eq)]
pub struct PublicKey([u8; 32]);

/// An Ed25519 signature (64 bytes: `R || S`).
#[derive(Clone, PartialEq, Eq)]
pub struct Signature([u8; 64]);

/// An Ed25519 keypair.
///
/// Holds the expanded secret (the scalar and nonce prefix) and the public key,
/// all derived once at construction, so [`Keypair::sign`] needs neither to
/// re-hash the seed nor to recompute the public key.
#[derive(Clone)]
pub struct Keypair {
    secret: SecretKey,
    public: PublicKey,
    // expanded secret, cached for fast repeated signing
    scalar: Scalar,
    prefix: [u8; 32],
}

impl SecretKey {
    pub fn from_bytes(seed: [u8; 32]) -> Self {
        SecretKey(seed)
    }
    pub fn to_bytes(&self) -> [u8; 32] {
        self.0
    }
    /// Derive the matching public key.
    pub fn public_key(&self) -> PublicKey {
        self.public_key_with(sha512)
    }

    /// Derive the matching public key, with the given hash function.
    ///
    /// This is the generic version of [`SecretKey::public_key`]; see the module
    /// documentation about alternative hash functions.
    pub fn public_key_with<H>(&self, hash: H) -> PublicKey
    where
        H: Fn(&[&[u8]]) -> [u8; HASH_LENGTH],
    {
        PublicKey(public_from_seed(&hash, &self.0))
    }

    /// Sign a message.
    pub fn sign(&self, message: &[u8]) -> Signature {
        self.sign_with(sha512, message)
    }

    /// Sign a message, with the given hash function.
    ///
    /// This is the generic version of [`SecretKey::sign`]; see the module
    /// documentation about alternative hash functions.
    pub fn sign_with<H>(&self, hash: H, message: &[u8]) -> Signature
    where
        H: Fn(&[&[u8]]) -> [u8; HASH_LENGTH],
    {
        Signature(sign(&hash, &self.0, message))
    }
}

impl PublicKey {
    pub fn from_bytes(bytes: [u8; 32]) -> Self {
        PublicKey(bytes)
    }
    pub fn to_bytes(&self) -> [u8; 32] {
        self.0
    }

    /// Verify a signature over `message`.
    pub fn verify(&self, message: &[u8], signature: &Signature) -> bool {
        self.verify_with(sha512, message, signature)
    }

    /// Verify a signature over `message`, with the given hash function.
    ///
    /// This is the generic version of [`PublicKey::verify`]; see the module
    /// documentation about alternative hash functions.
    pub fn verify_with<H>(&self, hash: H, message: &[u8], signature: &Signature) -> bool
    where
        H: Fn(&[&[u8]]) -> [u8; HASH_LENGTH],
    {
        verify(&hash, &self.0, message, &signature.0)
    }

    /// Do the part of verification that depends on the key alone, once, for
    /// verifying many signatures under this key.
    ///
    /// `None` if the key is not a valid point encoding
    pub fn precompute(&self) -> Option<PrecomputedPublicKey> {
        PrecomputedPublicKey::new(self)
    }
}

/// An Ed25519 public key with the per-key half of verification already done.
///
/// [`PublicKey::verify`] starts every signature by decompressing the 32 bytes
/// of the key into a curve point (field square root + builds a small table
/// of the multiples of that point for the double-scalar
/// multiplication).
///
/// A `PrecomputedPublicKey` pays for this once and use wider table, and
/// allow faster verification of multiple ed25519 signatures for different message
/// for the *same* public key.
///
/// This is a *precomputation*, not batch verification: signatures are still
/// verified one by one, with exactly the same accept/reject decision as
/// [`PublicKey::verify`] makes.
///
/// ```
/// # use eccoxide::protocol::ed25519::Keypair;
/// let keypair = Keypair::from_seed([42u8; 32]);
/// let messages: [&[u8]; 2] = [b"first", b"second"];
/// let signatures = messages.map(|m| keypair.sign(m));
///
/// let public = keypair.public().precompute().expect("a key we just built");
/// for (message, signature) in messages.iter().zip(signatures.iter()) {
///     assert!(public.verify(message, signature));
/// }
/// ```
#[derive(Clone)]
pub struct PrecomputedPublicKey {
    /// the key as it sits on the wire: `k = H(R || A || M)` hashes the
    /// encoding, not the point, so the bytes are needed as well as the table
    encoded: [u8; 32],
    /// the multiples of `-A`, the form the verification equation reads
    neg_multiples: PrecomputedPoint,
}

impl PrecomputedPublicKey {
    /// Precompute the verification tables of `public`, or `None` if it is not a
    /// valid point encoding.
    pub fn new(public: &PublicKey) -> Option<Self> {
        let a_point = decode_point(&public.0)?;
        Some(PrecomputedPublicKey {
            encoded: public.0,
            neg_multiples: PrecomputedPoint::new(&-&a_point),
        })
    }

    /// The public key the precomputation was built from.
    pub fn public_key(&self) -> PublicKey {
        PublicKey(self.encoded)
    }

    pub fn to_bytes(&self) -> [u8; 32] {
        self.encoded
    }

    /// Verify a signature over `message`, accepting exactly what
    /// [`PublicKey::verify`] accepts.
    pub fn verify(&self, message: &[u8], signature: &Signature) -> bool {
        self.verify_with(sha512, message, signature)
    }

    /// Verify a signature over `message` with the given hash function,
    /// accepting exactly what [`PublicKey::verify_with`] accepts.
    ///
    /// The precomputation itself is hash agnostic — it is the point
    /// decompression and the multiples of `-A`, neither of which hashes — so
    /// the same [`PrecomputedPublicKey`] serves any hash function.
    pub fn verify_with<H>(&self, hash: H, message: &[u8], signature: &Signature) -> bool
    where
        H: Fn(&[&[u8]]) -> [u8; HASH_LENGTH],
    {
        verify_equation(&hash, &self.encoded, message, &signature.0, |s, k| {
            self.neg_multiples.double_scalar_mul_base_vartime(s, k)
        })
    }
}

impl Signature {
    pub fn from_bytes(bytes: [u8; 64]) -> Self {
        Signature(bytes)
    }
    pub fn to_bytes(&self) -> [u8; 64] {
        self.0
    }
}

impl Keypair {
    /// Build a keypair from a 32-byte seed.
    pub fn from_seed(seed: [u8; 32]) -> Self {
        Self::from_seed_with(sha512, seed)
    }

    /// Build a keypair from a 32-byte seed, with the given hash function.
    ///
    /// This is the generic version of [`Keypair::from_seed`]; see the module
    /// documentation about alternative hash functions.
    ///
    /// The cached expanded secret is the one the hash function derives, so
    /// [`Keypair::sign_with`] must be handed the same hash function: signing
    /// such a keypair with another one mixes the two schemes.
    pub fn from_seed_with<H>(hash: H, seed: [u8; 32]) -> Self
    where
        H: Fn(&[&[u8]]) -> [u8; HASH_LENGTH],
    {
        // expand the seed once, reusing it for both the public key and the
        // cached signing material
        let (scalar, prefix) = expand_secret(&hash, &seed);
        let public = PublicKey(encode_point(&Point::mul_base(&scalar)));
        Keypair {
            secret: SecretKey::from_bytes(seed),
            public,
            scalar,
            prefix,
        }
    }

    pub fn public(&self) -> &PublicKey {
        &self.public
    }
    pub fn secret(&self) -> &SecretKey {
        &self.secret
    }

    /// Sign a message.
    ///
    /// Uses the cached expanded secret and public key, so this performs a single
    /// fixed-base scalar multiplication (for the nonce point R) instead of the
    /// two that signing from a bare seed would need.
    pub fn sign(&self, message: &[u8]) -> Signature {
        self.sign_with(sha512, message)
    }

    /// Sign a message, with the given hash function.
    ///
    /// This is the generic version of [`Keypair::sign`]; see the module
    /// documentation about alternative hash functions. The hash function must
    /// be the one [`Keypair::from_seed_with`] expanded the seed with.
    pub fn sign_with<H>(&self, hash: H, message: &[u8]) -> Signature
    where
        H: Fn(&[&[u8]]) -> [u8; HASH_LENGTH],
    {
        Signature(sign_with_public(
            &hash,
            &self.scalar,
            &self.prefix,
            &self.public.0,
            message,
        ))
    }
}

/// Ed25519 using BLAKE2b-512 as its hash function
///
/// Every type here mirrors the identically named type of the parent module,
/// with [`blake2b512`] substituted for [`sha512`]. This is *not* the standard
/// scheme of RFC 8032: keys and signatures produced here are only compatible
/// with this variant, which is why they are types of their own — the two
/// cannot be mixed up.
///
/// ```
/// use eccoxide::protocol::ed25519::blake2b::Keypair;
///
/// let keypair = Keypair::from_seed([42u8; 32]);
/// let signature = keypair.sign(b"message");
/// assert!(keypair.public().verify(b"message", &signature));
/// ```
#[cfg(feature = "ed25519-blake2")]
pub mod blake2b {
    use super::blake2b512;

    /// An Ed25519 signature (64 bytes: `R || S`)
    ///
    /// The signature format does not depend on the hash function, so this is
    /// the [`Signature`] of the parent module.
    pub use super::Signature;

    /// An Ed25519-BLAKE2b secret key: the 32-byte seed
    ///
    /// See [`SecretKey`](super::SecretKey).
    #[derive(Clone)]
    pub struct SecretKey(super::SecretKey);

    /// An Ed25519-BLAKE2b public key (32-byte compressed point)
    ///
    /// See [`PublicKey`](super::PublicKey).
    #[derive(Clone, PartialEq, Eq)]
    pub struct PublicKey(super::PublicKey);

    /// An Ed25519-BLAKE2b keypair
    ///
    /// See [`Keypair`](super::Keypair).
    #[derive(Clone)]
    pub struct Keypair(super::Keypair);

    /// An Ed25519-BLAKE2b public key with the per-key half of verification
    /// already done
    ///
    /// See [`PrecomputedPublicKey`](super::PrecomputedPublicKey).
    #[derive(Clone)]
    pub struct PrecomputedPublicKey(super::PrecomputedPublicKey);

    impl SecretKey {
        pub fn from_bytes(seed: [u8; 32]) -> Self {
            SecretKey(super::SecretKey::from_bytes(seed))
        }
        pub fn to_bytes(&self) -> [u8; 32] {
            self.0.to_bytes()
        }
        /// Derive the matching public key.
        pub fn public_key(&self) -> PublicKey {
            PublicKey(self.0.public_key_with(blake2b512))
        }
        /// Sign a message.
        pub fn sign(&self, message: &[u8]) -> Signature {
            self.0.sign_with(blake2b512, message)
        }
    }

    impl PublicKey {
        pub fn from_bytes(bytes: [u8; 32]) -> Self {
            PublicKey(super::PublicKey::from_bytes(bytes))
        }

        pub fn to_bytes(&self) -> [u8; 32] {
            self.0.to_bytes()
        }

        /// Verify a signature over `message`.
        pub fn verify(&self, message: &[u8], signature: &Signature) -> bool {
            self.0.verify_with(blake2b512, message, signature)
        }

        /// Do the part of verification that depends on the key alone, once, for
        /// verifying many signatures under this key.
        ///
        /// `None` if the key is not a valid point encoding
        pub fn precompute(&self) -> Option<PrecomputedPublicKey> {
            PrecomputedPublicKey::new(self)
        }
    }

    impl PrecomputedPublicKey {
        /// Precompute the verification tables of `public`, or `None` if it is
        /// not a valid point encoding.
        pub fn new(public: &PublicKey) -> Option<Self> {
            super::PrecomputedPublicKey::new(&public.0).map(PrecomputedPublicKey)
        }

        /// The public key the precomputation was built from.
        pub fn public_key(&self) -> PublicKey {
            PublicKey(self.0.public_key())
        }

        pub fn to_bytes(&self) -> [u8; 32] {
            self.0.to_bytes()
        }

        /// Verify a signature over `message`, accepting exactly what
        /// [`PublicKey::verify`] accepts.
        pub fn verify(&self, message: &[u8], signature: &Signature) -> bool {
            self.0.verify_with(blake2b512, message, signature)
        }
    }

    impl Keypair {
        /// Build a keypair from a 32-byte seed.
        pub fn from_seed(seed: [u8; 32]) -> Self {
            Keypair(super::Keypair::from_seed_with(blake2b512, seed))
        }

        /// The public key of the keypair.
        ///
        /// Unlike [`Keypair::public`](super::Keypair::public) this returns by
        /// value: the key the parent keypair holds is the SHA-512 type, and the
        /// two variants are distinct types.
        pub fn public(&self) -> PublicKey {
            PublicKey(self.0.public().clone())
        }

        /// The secret key of the keypair, returned by value for the same reason
        /// as [`Keypair::public`].
        pub fn secret(&self) -> SecretKey {
            SecretKey(self.0.secret().clone())
        }

        /// Sign a message.
        pub fn sign(&self, message: &[u8]) -> Signature {
            self.0.sign_with(blake2b512, message)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hex<const N: usize>(s: &str) -> [u8; N] {
        let mut out = [0u8; N];
        for (i, b) in out.iter_mut().enumerate() {
            *b = u8::from_str_radix(&s[i * 2..i * 2 + 2], 16).unwrap();
        }
        out
    }

    struct Vector {
        seed: &'static str,
        public: &'static str,
        message: &'static str,
        signature: &'static str,
    }

    // RFC 8032, section 7.1 (TEST 1, 2, 3). Signatures cross-checked against an
    // independent RFC-8032 reference implementation reproducing the published
    // public keys.
    const VECTORS: &[Vector] = &[
        Vector {
            seed: "9d61b19deffd5a60ba844af492ec2cc44449c5697b326919703bac031cae7f60",
            public: "d75a980182b10ab7d54bfed3c964073a0ee172f3daa62325af021a68f707511a",
            message: "",
            signature: "e5564300c360ac729086e2cc806e828a84877f1eb8e5d974d873e065224901555fb8821590a33bacc61e39701cf9b46bd25bf5f0595bbe24655141438e7a100b",
        },
        Vector {
            seed: "4ccd089b28ff96da9db6c346ec114e0f5b8a319f35aba624da8cf6ed4fb8a6fb",
            public: "3d4017c3e843895a92b70aa74d1b7ebc9c982ccf2ec4968cc0cd55f12af4660c",
            message: "72",
            signature: "92a009a9f0d4cab8720e820b5f642540a2b27b5416503f8fb3762223ebdb69da085ac1e43e15996e458f3613d0f11d8c387b2eaeb4302aeeb00d291612bb0c00",
        },
        Vector {
            seed: "c5aa8df43f9f837bedb7442f31dcb7b166d38535076f094b85ce3a2e0b4458f7",
            public: "fc51cd8e6218a1a38da47ed00230f0580816ed13ba3303ac5deb911548908025",
            message: "af82",
            signature: "6291d657deec24024827e69c3abe01a30ce548a284743a445e3680d7db5ac3ac18ff9b538d16f290ae67f760984dc6594a7c15e9716ed28dc027beceea1ec40a",
        },
    ];

    fn hex_vec(s: &str) -> Vec<u8> {
        (0..s.len() / 2)
            .map(|i| u8::from_str_radix(&s[i * 2..i * 2 + 2], 16).unwrap())
            .collect()
    }

    #[test]
    fn rfc8032_vectors() {
        for v in VECTORS {
            let seed: [u8; 32] = hex(v.seed);
            let expected_pub: [u8; 32] = hex(v.public);
            let expected_sig: [u8; 64] = hex(v.signature);
            let message = hex_vec(v.message);

            let sk = SecretKey::from_bytes(seed);
            let pk = sk.public_key();
            assert_eq!(pk.to_bytes(), expected_pub, "public key mismatch");

            let sig = sk.sign(&message);
            assert_eq!(sig.to_bytes(), expected_sig, "signature mismatch");

            assert!(pk.verify(&message, &sig), "valid signature rejected");
        }
    }

    #[test]
    fn keypair_sign_matches_secretkey_and_rfc() {
        // The cached-public-key Keypair path must produce byte-identical
        // signatures to the from-seed SecretKey path (Ed25519 is deterministic),
        // and both must match the RFC 8032 vectors.
        for v in VECTORS {
            let seed: [u8; 32] = hex(v.seed);
            let expected_sig: [u8; 64] = hex(v.signature);
            let message = hex_vec(v.message);

            let kp = Keypair::from_seed(seed);
            let sk = SecretKey::from_bytes(seed);

            assert_eq!(kp.public().to_bytes(), hex::<32>(v.public));
            let kp_sig = kp.sign(&message);
            assert_eq!(kp_sig.to_bytes(), sk.sign(&message).to_bytes());
            assert_eq!(kp_sig.to_bytes(), expected_sig, "keypair sig != RFC");
        }
    }

    /// Decoding written the straightforward way, recovering x and testing it
    /// for zero, as the reference the y-based rejection has to match.
    fn decode_point_reference(bytes: &[u8; 32]) -> Option<Point> {
        let x_sign_bit = bytes[31] >> 7;
        let mut le = *bytes;
        le[31] &= 0x7f;
        let y = FieldElement::from_bytes_le(&le)?;
        let want = if x_sign_bit == 1 {
            Sign::Negative
        } else {
            Sign::Positive
        };
        let p = Point::decompress(&y, want)?;
        let (x, _) = p.to_affine();
        if x_sign_bit == 1 && x.is_zero() {
            return None;
        }
        Some(p)
    }

    #[test]
    fn decode_point_matches_reference() {
        let mut cases: Vec<[u8; 32]> = Vec::new();

        // the two curve points with x = 0: y = 1 and y = p - 1
        let mut y_one = [0u8; 32];
        y_one[0] = 1;
        let mut y_minus_one = [0xffu8; 32];
        y_minus_one[0] = 0xec;
        y_minus_one[31] = 0x7f;
        for base in [y_one, y_minus_one] {
            cases.push(base);
            let mut sign_set = base;
            sign_set[31] |= 0x80; // the non-canonical encodings
            cases.push(sign_set);
        }

        // genuine points, with the sign bit left alone, forced and cleared
        for k in 1..20u64 {
            let encoded = encode_point(&Point::mul_base(&Scalar::from_u64(k)));
            cases.push(encoded);
            let mut sign_set = encoded;
            sign_set[31] |= 0x80;
            cases.push(sign_set);
            let mut sign_clear = encoded;
            sign_clear[31] &= 0x7f;
            cases.push(sign_clear);
        }

        // pseudo-random strings: mostly off-curve or non-canonical y
        let mut h = sha512(&[b"ed25519 decode cases"]);
        for i in 0..200u8 {
            let mut bytes = [0u8; 32];
            bytes.copy_from_slice(&h[..32]);
            cases.push(bytes);
            h = sha512(&[&h[..], &[i]]);
        }

        let mut accepted = 0;
        for bytes in cases {
            let decoded = decode_point(&bytes);
            accepted += decoded.is_some() as usize;
            assert_eq!(
                decoded,
                decode_point_reference(&bytes),
                "decoding disagrees for {:02x?}",
                bytes
            );
        }
        // the cases do cover valid encodings, not only rejections
        assert!(accepted > 20, "only {} accepted", accepted);
    }

    #[test]
    fn tampered_signature_and_message_rejected() {
        let kp = Keypair::from_seed(hex(VECTORS[1].seed));
        let msg = b"attack at dawn";
        let sig = kp.sign(msg);
        assert!(kp.public().verify(msg, &sig));

        // wrong message
        assert!(!kp.public().verify(b"attack at dusk", &sig));

        // flipped bit in the signature
        let mut bad = sig.to_bytes();
        bad[10] ^= 1;
        assert!(!kp.public().verify(msg, &Signature::from_bytes(bad)));

        // wrong public key
        let other = Keypair::from_seed(hex(VECTORS[2].seed));
        assert!(!other.public().verify(msg, &sig));
    }

    #[test]
    fn precomputed_verify_matches_rfc_vectors() {
        for v in VECTORS {
            let pk = PublicKey::from_bytes(hex(v.public));
            let sig = Signature::from_bytes(hex(v.signature));
            let message = hex_vec(v.message);
            let precomputed = pk.precompute().expect("RFC public key decodes");

            assert_eq!(precomputed.to_bytes(), pk.to_bytes());
            assert!(precomputed.public_key() == pk);
            assert!(
                precomputed.verify(&message, &sig),
                "valid signature rejected"
            );
        }
    }

    #[test]
    fn precomputed_verify_agrees_with_plain_verify() {
        // the precomputed path must reach the same decision as the plain one on
        // every case, accepted or rejected: the tables are an optimisation, not
        // a different check
        let keypairs: Vec<Keypair> = VECTORS
            .iter()
            .map(|v| Keypair::from_seed(hex(v.seed)))
            .chain(core::iter::once(Keypair::from_seed([3u8; 32])))
            .collect();

        for kp in keypairs.iter() {
            let precomputed = kp.public().precompute().expect("derived key decodes");
            for len in [0usize, 1, 32, 33, 150] {
                let message: Vec<u8> = (0..len).map(|i| (i * 7 + 5) as u8).collect();
                let sig = kp.sign(&message);

                let mut cases = alloc::vec![
                    // the genuine signature
                    (message.clone(), sig.to_bytes()),
                    // a different message under the same signature
                    (alloc::vec![0xaa; len + 1], sig.to_bytes()),
                ];
                // every kind of damage to the 64 signature bytes: R, S, the
                // sign bit of R, and an S past the group order
                for byte in [0usize, 31, 32, 63] {
                    let mut bad = sig.to_bytes();
                    bad[byte] ^= 0x80;
                    cases.push((message.clone(), bad));
                }
                let mut s_max = sig.to_bytes();
                s_max[32..].copy_from_slice(&[0xff; 32]); // S >= l, non-canonical
                cases.push((message.clone(), s_max));

                for (m, sig_bytes) in cases {
                    let signature = Signature::from_bytes(sig_bytes);
                    assert_eq!(
                        precomputed.verify(&m, &signature),
                        kp.public().verify(&m, &signature),
                        "verdicts differ for a {}-byte message",
                        m.len()
                    );
                }

                // and under the wrong key
                for other in keypairs.iter() {
                    let other_precomputed =
                        other.public().precompute().expect("derived key decodes");
                    assert_eq!(
                        other_precomputed.verify(&message, &sig),
                        other.public().verify(&message, &sig),
                    );
                }
            }
        }
    }

    #[test]
    fn precompute_rejects_what_verify_rejects() {
        // a key that cannot be decoded has no precomputation, which is exactly
        // the case where verifying with it always fails
        let kp = Keypair::from_seed([9u8; 32]);
        let message = b"precomputation";
        let sig = kp.sign(message);

        let mut bad_keys: Vec<[u8; 32]> = alloc::vec![
            // x = 0 with the sign bit set: the non-canonical encodings
            {
                let mut y_one = [0u8; 32];
                y_one[0] = 1;
                y_one[31] = 0x80;
                y_one
            },
        ];
        // strings that are almost always off-curve
        let mut h = sha512(&[b"ed25519 precompute rejection"]);
        for i in 0..32u8 {
            let mut bytes = [0u8; 32];
            bytes.copy_from_slice(&h[..32]);
            bad_keys.push(bytes);
            h = sha512(&[&h[..], &[i]]);
        }

        let mut rejected = 0;
        for bytes in bad_keys {
            let pk = PublicKey::from_bytes(bytes);
            match pk.precompute() {
                None => {
                    rejected += 1;
                    assert!(!pk.verify(message, &sig), "undecodable key still verifies");
                }
                Some(precomputed) => {
                    // a decodable key: it is simply not the signer
                    assert_eq!(precomputed.verify(message, &sig), pk.verify(message, &sig));
                }
            }
        }
        assert!(rejected > 10, "only {} keys rejected", rejected);
    }

    #[test]
    fn roundtrip_various_messages() {
        let kp = Keypair::from_seed([7u8; 32]);
        for len in [0usize, 1, 31, 32, 33, 64, 200] {
            let msg: Vec<u8> = (0..len).map(|i| (i * 3 + 1) as u8).collect();
            let sig = kp.sign(&msg);
            assert!(kp.public().verify(&msg, &sig), "len={}", len);
        }
    }

    #[test]
    fn generic_hash_with_sha512_matches_plain() {
        // handed the standard hash, the parametrized interface must reproduce
        // the plain API exactly, RFC 8032 vectors included
        for v in VECTORS {
            let seed: [u8; 32] = hex(v.seed);
            let expected_sig: [u8; 64] = hex(v.signature);
            let message = hex_vec(v.message);

            let sk = SecretKey::from_bytes(seed);
            assert_eq!(sk.public_key_with(sha512).to_bytes(), hex::<32>(v.public));
            assert_eq!(sk.sign_with(sha512, &message).to_bytes(), expected_sig);

            let kp = Keypair::from_seed_with(sha512, seed);
            assert_eq!(kp.public().to_bytes(), hex::<32>(v.public));
            assert_eq!(kp.sign_with(sha512, &message).to_bytes(), expected_sig);

            let sig = Signature::from_bytes(expected_sig);
            let pk = PublicKey::from_bytes(hex(v.public));
            assert!(pk.verify_with(sha512, &message, &sig));
            assert!(
                pk.precompute()
                    .expect("RFC public key decodes")
                    .verify_with(sha512, &message, &sig)
            );
        }
    }

    /// A hash function of our own: SHA-512 with a domain separation prefix.
    /// Nothing standard — the point is that an arbitrary
    /// `Fn(&[&[u8]]) -> [u8; HASH_LENGTH]` is usable and yields a scheme that
    /// is consistent with itself and with nothing else.
    fn prefixed_sha512(slices: &[&[u8]]) -> [u8; HASH_LENGTH] {
        let mut parts: Vec<&[u8]> = alloc::vec![b"eccoxide test hash"];
        parts.extend_from_slice(slices);
        sha512(&parts)
    }

    #[test]
    fn arbitrary_hash_is_self_consistent() {
        let seed = [11u8; 32];
        let message = b"custom hash";

        let kp = Keypair::from_seed_with(prefixed_sha512, seed);
        let sig = kp.sign_with(prefixed_sha512, message);

        // the same scheme reached from the bare secret key
        let sk = SecretKey::from_bytes(seed);
        assert_eq!(
            sk.public_key_with(prefixed_sha512).to_bytes(),
            kp.public().to_bytes()
        );
        assert_eq!(
            sk.sign_with(prefixed_sha512, message).to_bytes(),
            sig.to_bytes()
        );

        assert!(kp.public().verify_with(prefixed_sha512, message, &sig));
        assert!(
            kp.public()
                .precompute()
                .expect("derived key decodes")
                .verify_with(prefixed_sha512, message, &sig)
        );

        // and a scheme of its own: the standard one neither derives the same
        // key nor accepts the signature
        let standard = Keypair::from_seed(seed);
        assert_ne!(kp.public().to_bytes(), standard.public().to_bytes());
        assert!(!kp.public().verify(message, &sig));
    }

    /// The BLAKE2b-512 variant, cross-checked against the reference
    /// implementation of `cryptoxide::ed25519::blake2b`.
    #[cfg(feature = "ed25519-blake2")]
    mod blake2b_variant {
        use super::super::{Keypair, PublicKey, Signature, blake2b, blake2b512};
        use cryptoxide::ed25519::blake2b as reference;

        const SEEDS: [[u8; 32]; 3] = [[0u8; 32], [1u8; 32], [0xabu8; 32]];
        const MESSAGES: [&[u8]; 4] = [b"", b"a", b"the quick brown fox", &[0xffu8; 100]];

        #[test]
        fn matches_cryptoxide() {
            for seed in SEEDS {
                let (ref_keypair, ref_public) = reference::keypair(&seed);
                let kp = blake2b::Keypair::from_seed(seed);
                let sk = blake2b::SecretKey::from_bytes(seed);

                assert_eq!(kp.public().to_bytes(), ref_public);
                assert_eq!(sk.public_key().to_bytes(), ref_public);
                assert_eq!(kp.secret().to_bytes(), seed);

                for message in MESSAGES {
                    let sig = kp.sign(message);
                    assert_eq!(sig.to_bytes(), reference::signature(message, &ref_keypair));
                    assert_eq!(sk.sign(message).to_bytes(), sig.to_bytes());
                    assert!(reference::verify(message, &ref_public, &sig.to_bytes()));
                    assert!(kp.public().verify(message, &sig));
                }
            }
        }

        #[test]
        fn is_the_parametrized_scheme_on_blake2b512() {
            // the module is exactly the generic interface pre-applied
            let seed = [3u8; 32];
            let module = blake2b::Keypair::from_seed(seed);
            let generic = Keypair::from_seed_with(blake2b512, seed);
            assert_eq!(module.public().to_bytes(), generic.public().to_bytes());
            for message in MESSAGES {
                assert_eq!(
                    module.sign(message).to_bytes(),
                    generic.sign_with(blake2b512, message).to_bytes()
                );
            }
        }

        #[test]
        fn precomputed_verify_agrees_with_plain() {
            let kp = blake2b::Keypair::from_seed([5u8; 32]);
            let public = kp.public();
            let precomputed = public.precompute().expect("derived key decodes");
            assert_eq!(precomputed.to_bytes(), public.to_bytes());
            assert!(precomputed.public_key() == public);

            for message in MESSAGES {
                let sig = kp.sign(message);
                assert!(precomputed.verify(message, &sig));
                assert!(public.verify(message, &sig));

                // damaged R, damaged S, and a non-canonical S
                let mut cases = alloc::vec![];
                for byte in [0usize, 31, 32, 63] {
                    let mut bad = sig.to_bytes();
                    bad[byte] ^= 0x80;
                    cases.push(bad);
                }
                let mut s_max = sig.to_bytes();
                s_max[32..].copy_from_slice(&[0xff; 32]);
                cases.push(s_max);

                for bad in cases {
                    let bad = Signature::from_bytes(bad);
                    assert!(!public.verify(message, &bad));
                    assert_eq!(
                        precomputed.verify(message, &bad),
                        public.verify(message, &bad)
                    );
                }
            }
        }

        #[test]
        fn not_interchangeable_with_sha512() {
            let seed = [7u8; 32];
            let message: &[u8] = b"hash binding";
            let variant = blake2b::Keypair::from_seed(seed);
            let standard = Keypair::from_seed(seed);
            assert_ne!(variant.public().to_bytes(), standard.public().to_bytes());

            // the standard scheme rejects a BLAKE2b signature, under its own
            // key and under the BLAKE2b one
            let variant_sig = variant.sign(message);
            assert!(!standard.public().verify(message, &variant_sig));
            assert!(
                !PublicKey::from_bytes(variant.public().to_bytes()).verify(message, &variant_sig)
            );

            // and the variant rejects a standard signature
            let standard_sig = standard.sign(message);
            assert!(!variant.public().verify(message, &standard_sig));
            let standard_as_variant = blake2b::PublicKey::from_bytes(standard.public().to_bytes());
            assert!(!standard_as_variant.verify(message, &standard_sig));
        }
    }
}
