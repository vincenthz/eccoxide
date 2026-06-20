//! Ed25519 signature scheme (RFC 8032, PureEdDSA).
//!
//! Keys, signatures and encoded points use the RFC 8032 little-endian wire
//! format. The internal SHA-512 is provided by the `cryptoxide` crate.

use crate::curve::curve25519::{FieldElement, Point, Scalar};
use crate::curve::field::Sign;
use cryptoxide::hashing::sha2::Sha512;

/// SHA-512 of the concatenation of `parts`.
fn sha512(parts: &[&[u8]]) -> [u8; 64] {
    let mut h = Sha512::new();
    for p in parts {
        h.update_mut(p);
    }
    h.finalize()
}

/// Reduce a 64-byte SHA-512 output, interpreted little-endian, modulo the group
/// order `l`.
fn reduce_wide_le(hash: &[u8; 64]) -> Scalar {
    let mut be = *hash;
    be.reverse(); // the wide reducer is big-endian
    Scalar::init_from_wide_bytes(be)
}

/// Encode a point to its 32-byte RFC 8032 representation: little-endian y with
/// the low bit of x stored in the most significant bit.
fn encode_point(p: &Point) -> [u8; 32] {
    let (x, y) = p.to_affine();
    let mut out = y.to_bytes(); // big-endian
    out.reverse(); // little-endian
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
    le.reverse(); // big-endian
    // y must be a canonical field element (< p)
    let y = FieldElement::from_bytes(&le)?;
    let want = if x_sign_bit == 1 {
        Sign::Negative
    } else {
        Sign::Positive
    };
    let p = Point::decompress(&y, want)?;
    // reject the (x = 0, sign = 1) non-canonical encoding
    let (x, _) = p.to_affine();
    if x_sign_bit == 1 && x.is_zero() {
        return None;
    }
    Some(p)
}

/// Expand a 32-byte seed into the secret scalar (mod l) and the nonce prefix.
fn expand_secret(seed: &[u8; 32]) -> (Scalar, [u8; 32]) {
    let h = sha512(&[&seed[..]]);

    let mut a_le = [0u8; 32];
    a_le.copy_from_slice(&h[..32]);
    a_le[0] &= 248;
    a_le[31] &= 127;
    a_le[31] |= 64;

    // reduce the clamped little-endian scalar modulo l
    let mut a_be = a_le;
    a_be.reverse();
    let mut wide = [0u8; 64];
    wide[32..].copy_from_slice(&a_be);
    let a = Scalar::init_from_wide_bytes(wide);

    let mut prefix = [0u8; 32];
    prefix.copy_from_slice(&h[32..]);
    (a, prefix)
}

fn public_from_seed(seed: &[u8; 32]) -> [u8; 32] {
    let (a, _) = expand_secret(seed);
    encode_point(&Point::GENERATOR.scale(&a))
}

fn sign(seed: &[u8; 32], message: &[u8]) -> [u8; 64] {
    let (a, prefix) = expand_secret(seed);
    let public = encode_point(&Point::GENERATOR.scale(&a));

    // r = H(prefix || M) mod l ; R = [r]B
    let r = reduce_wide_le(&sha512(&[&prefix[..], message]));
    let r_encoded = encode_point(&Point::GENERATOR.scale(&r));

    // k = H(R || A || M) mod l
    let k = reduce_wide_le(&sha512(&[&r_encoded[..], &public[..], message]));

    // S = (r + k·a) mod l
    let s = &r + &(&k * &a);
    let mut s_le = s.to_bytes(); // big-endian
    s_le.reverse(); // little-endian

    let mut sig = [0u8; 64];
    sig[..32].copy_from_slice(&r_encoded);
    sig[32..].copy_from_slice(&s_le);
    sig
}

fn verify(public: &[u8; 32], message: &[u8], sig: &[u8; 64]) -> bool {
    let a_point = match decode_point(public) {
        Some(p) => p,
        None => return false,
    };
    let mut r_encoded = [0u8; 32];
    r_encoded.copy_from_slice(&sig[..32]);
    let r_point = match decode_point(&r_encoded) {
        Some(p) => p,
        None => return false,
    };

    // S must be canonical (< l)
    let mut s_be = [0u8; 32];
    s_be.copy_from_slice(&sig[32..]);
    s_be.reverse();
    let s = match Scalar::from_bytes(&s_be) {
        Some(s) => s,
        None => return false,
    };

    // k = H(R || A || M) mod l
    let k = reduce_wide_le(&sha512(&[&r_encoded[..], &public[..], message]));

    // accept iff [S]B == R + [k]A
    let lhs = Point::GENERATOR.scale(&s);
    let rhs = &r_point + &a_point.scale(&k);
    lhs == rhs
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
#[derive(Clone)]
pub struct Keypair {
    secret: SecretKey,
    public: PublicKey,
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
        PublicKey(public_from_seed(&self.0))
    }
    /// Sign a message.
    pub fn sign(&self, message: &[u8]) -> Signature {
        Signature(sign(&self.0, message))
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
        verify(&self.0, message, &signature.0)
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
        let secret = SecretKey::from_bytes(seed);
        let public = secret.public_key();
        Keypair { secret, public }
    }
    pub fn public(&self) -> &PublicKey {
        &self.public
    }
    pub fn secret(&self) -> &SecretKey {
        &self.secret
    }
    pub fn sign(&self, message: &[u8]) -> Signature {
        self.secret.sign(message)
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
    fn roundtrip_various_messages() {
        let kp = Keypair::from_seed([7u8; 32]);
        for len in [0usize, 1, 31, 32, 33, 64, 200] {
            let msg: Vec<u8> = (0..len).map(|i| (i * 3 + 1) as u8).collect();
            let sig = kp.sign(&msg);
            assert!(kp.public().verify(&msg, &sig), "len={}", len);
        }
    }
}
