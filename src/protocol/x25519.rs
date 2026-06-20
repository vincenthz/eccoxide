//! X25519 Diffie-Hellman key agreement (RFC 7748).
//!
//! This is the Montgomery-ladder ECDH function on curve25519. Scalars and
//! u-coordinates are exchanged as 32-byte little-endian arrays, as specified by
//! the RFC.

use crate::curve::curve25519::{FieldElement, MontgomeryPoint};

/// The standard base point u-coordinate (the integer 9), little-endian.
pub const BASEPOINT: [u8; 32] = [
    9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

/// `decodeScalar25519`: clamp the scalar as required by X25519.
fn clamp(mut k: [u8; 32]) -> [u8; 32] {
    k[0] &= 248;
    k[31] &= 127;
    k[31] |= 64;
    k
}

/// `decodeUCoordinate`: mask the unused top bit and read the u-coordinate as a
/// field element (the value is implicitly reduced modulo p by the arithmetic).
fn decode_u(u: &[u8; 32]) -> FieldElement {
    let mut le = *u;
    le[31] &= 0x7f;
    // curve25519 field elements are little-endian, matching the wire format
    FieldElement::from_bytes_unchecked_le(&le)
}

/// The raw X25519 function: multiply the u-coordinate point `u` by the (clamped)
/// scalar `scalar`, returning the resulting u-coordinate (little-endian).
///
/// A result of all-zero bytes indicates a low-order input point and should be
/// rejected by callers performing a Diffie-Hellman exchange.
pub fn x25519(scalar: &[u8; 32], u: &[u8; 32]) -> [u8; 32] {
    let mut k_be = clamp(*scalar);
    k_be.reverse(); // the ladder consumes a big-endian scalar

    let point = MontgomeryPoint::from_u(&decode_u(u));
    let result = point.scale_bytes(&k_be).u();

    // little-endian wire format, matching the native field byte order
    result.to_bytes_le()
}

/// X25519 against the standard base point: computes the public u-coordinate for
/// a secret `scalar`.
pub fn x25519_base(scalar: &[u8; 32]) -> [u8; 32] {
    x25519(scalar, &BASEPOINT)
}

/// An X25519 secret key (32 bytes of secret scalar material, clamped on use).
#[derive(Clone)]
pub struct SecretKey([u8; 32]);

/// An X25519 public key (a u-coordinate, little-endian).
#[derive(Clone, PartialEq, Eq)]
pub struct PublicKey([u8; 32]);

/// The shared secret produced by [`SecretKey::diffie_hellman`].
#[derive(Clone)]
pub struct SharedSecret([u8; 32]);

impl SecretKey {
    pub fn from_bytes(bytes: [u8; 32]) -> Self {
        SecretKey(bytes)
    }

    pub fn to_bytes(&self) -> [u8; 32] {
        self.0
    }

    /// Derive the matching public key.
    pub fn public_key(&self) -> PublicKey {
        PublicKey(x25519_base(&self.0))
    }

    /// Perform the Diffie-Hellman exchange with a peer's public key.
    pub fn diffie_hellman(&self, peer: &PublicKey) -> SharedSecret {
        SharedSecret(x25519(&self.0, &peer.0))
    }
}

impl PublicKey {
    pub fn from_bytes(bytes: [u8; 32]) -> Self {
        PublicKey(bytes)
    }
    pub fn to_bytes(&self) -> [u8; 32] {
        self.0
    }
}

impl SharedSecret {
    pub fn to_bytes(&self) -> [u8; 32] {
        self.0
    }

    /// True if the shared secret is the all-zero value indicating a low-order
    /// peer point. Callers should reject such exchanges.
    pub fn was_contributory(&self) -> bool {
        self.0 != [0u8; 32]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn h(s: &str) -> [u8; 32] {
        let mut out = [0u8; 32];
        for (i, b) in out.iter_mut().enumerate() {
            *b = u8::from_str_radix(&s[i * 2..i * 2 + 2], 16).unwrap();
        }
        out
    }

    // RFC 7748, section 5.2
    #[test]
    fn rfc7748_vector1() {
        let k = h("a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4");
        let u = h("e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c");
        let r = h("c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552");
        assert_eq!(x25519(&k, &u), r);
    }

    #[test]
    fn rfc7748_vector2() {
        let k = h("4b66e9d4d1b4673c5ad22691957d6af5c11b6421e0ea01d42ca4169e7918ba0d");
        let u = h("e5210f12786811d3f4b7959d0538ae2c31dbe7106fc03c3efc4cd549c715a493");
        let r = h("95cbde9476e8907d7aade45cb4b873f88b595a68799fa152e6f8f7647aac7957");
        assert_eq!(x25519(&k, &u), r);
    }

    // RFC 7748, section 5.2: one iteration of the base-point recurrence.
    #[test]
    fn rfc7748_iterated_once() {
        let k = h("0900000000000000000000000000000000000000000000000000000000000000");
        let r = h("422c8e7a6227d7bca1350b3e2bb7279f7897b87bb6854b783c60e80311ae3079");
        assert_eq!(x25519_base(&k), r);
    }

    #[test]
    fn diffie_hellman_agrees() {
        let a = SecretKey::from_bytes(h(
            "77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a",
        ));
        let b = SecretKey::from_bytes(h(
            "5dab087e624a8a4b79e17f8b83800ee66f3bb1292618b6fd1c2f8b27ff88e0eb",
        ));
        let ab = a.diffie_hellman(&b.public_key());
        let ba = b.diffie_hellman(&a.public_key());
        assert_eq!(ab.to_bytes(), ba.to_bytes());
        // RFC 7748 section 6.1 expected shared secret
        assert_eq!(
            ab.to_bytes(),
            h("4a5d9d5ba4ce2de1728e3bf480350f25e07e21c947d19e3376f09b3c1e161742")
        );
        assert!(ab.was_contributory());
    }
}
