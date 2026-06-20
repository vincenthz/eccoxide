//! X448 Diffie-Hellman key agreement (RFC 7748).
//!
//! This is the Montgomery-ladder ECDH function on curve448. Scalars and
//! u-coordinates are exchanged as 56-byte little-endian arrays, as specified by
//! the RFC.

use crate::curve::curve448::{FieldElement, MontgomeryPoint};

/// The standard base point u-coordinate (the integer 5), little-endian.
pub const BASEPOINT: [u8; 56] = [
    5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

/// `decodeScalar448`: clamp the scalar as required by X448.
fn clamp(mut k: [u8; 56]) -> [u8; 56] {
    k[0] &= 252;
    k[55] |= 128;
    k
}

/// `decodeUCoordinate`: read the u-coordinate as a field element. For X448 the
/// field is exactly 448 bits = 56 bytes, so (unlike X25519) no bit is masked.
fn decode_u(u: &[u8; 56]) -> FieldElement {
    // curve448 field elements are little-endian, matching the wire format
    FieldElement::from_bytes_unchecked_le(u)
}

/// The raw X448 function: multiply the u-coordinate point `u` by the (clamped)
/// scalar `scalar`, returning the resulting u-coordinate (little-endian).
///
/// A result of all-zero bytes indicates a low-order input point and should be
/// rejected by callers performing a Diffie-Hellman exchange.
pub fn x448(scalar: &[u8; 56], u: &[u8; 56]) -> [u8; 56] {
    let mut k_be = clamp(*scalar);
    k_be.reverse(); // the ladder consumes a big-endian scalar

    let point = MontgomeryPoint::from_u(&decode_u(u));
    let result = point.scale_bytes(&k_be).u();

    // little-endian wire format, matching the native field byte order
    result.to_bytes_le()
}

/// X448 against the standard base point: computes the public u-coordinate for
/// a secret `scalar`.
pub fn x448_base(scalar: &[u8; 56]) -> [u8; 56] {
    x448(scalar, &BASEPOINT)
}

/// An X448 secret key (56 bytes of secret scalar material, clamped on use).
#[derive(Clone)]
pub struct SecretKey([u8; 56]);

/// An X448 public key (a u-coordinate, little-endian).
#[derive(Clone, PartialEq, Eq)]
pub struct PublicKey([u8; 56]);

/// The shared secret produced by [`SecretKey::diffie_hellman`].
#[derive(Clone)]
pub struct SharedSecret([u8; 56]);

impl SecretKey {
    pub fn from_bytes(bytes: [u8; 56]) -> Self {
        SecretKey(bytes)
    }

    pub fn to_bytes(&self) -> [u8; 56] {
        self.0
    }

    /// Derive the matching public key.
    pub fn public_key(&self) -> PublicKey {
        PublicKey(x448_base(&self.0))
    }

    /// Perform the Diffie-Hellman exchange with a peer's public key.
    pub fn diffie_hellman(&self, peer: &PublicKey) -> SharedSecret {
        SharedSecret(x448(&self.0, &peer.0))
    }
}

impl PublicKey {
    pub fn from_bytes(bytes: [u8; 56]) -> Self {
        PublicKey(bytes)
    }
    pub fn to_bytes(&self) -> [u8; 56] {
        self.0
    }
}

impl SharedSecret {
    pub fn to_bytes(&self) -> [u8; 56] {
        self.0
    }

    /// True if the shared secret is the all-zero value indicating a low-order
    /// peer point. Callers should reject such exchanges.
    pub fn was_contributory(&self) -> bool {
        self.0 != [0u8; 56]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn h(s: &str) -> [u8; 56] {
        let mut out = [0u8; 56];
        for (i, b) in out.iter_mut().enumerate() {
            *b = u8::from_str_radix(&s[i * 2..i * 2 + 2], 16).unwrap();
        }
        out
    }

    // RFC 7748, section 5.2
    #[test]
    fn rfc7748_vector1() {
        let k = h("3d262fddf9ec8e88495266fea19a34d28882acef045104d0d1aae121700a779c984c24f8cdd78fbff44943eba368f54b29259a4f1c600ad3");
        let u = h("06fce640fa3487bfda5f6cf2d5263f8aad88334cbd07437f020f08f9814dc031ddbdc38c19c6da2583fa5429db94ada18aa7a7fb4ef8a086");
        let r = h("ce3e4ff95a60dc6697da1db1d85e6afbdf79b50a2412d7546d5f239fe14fbaadeb445fc66a01b0779d98223961111e21766282f73dd96b6f");
        assert_eq!(x448(&k, &u), r);
    }

    #[test]
    fn rfc7748_vector2() {
        let k = h("203d494428b8399352665ddca42f9de8fef600908e0d461cb021f8c538345dd77c3e4806e25f46d3315c44e0a5b4371282dd2c8d5be3095f");
        let u = h("0fbcc2f993cd56d3305b0b7d9e55d4c1a8fb5dbb52f8e9a1e9b6201b165d015894e56c4d3570bee52fe205e28a78b91cdfbde71ce8d157db");
        let r = h("884a02576239ff7a2f2f63b2db6a9ff37047ac13568e1e30fe63c4a7ad1b3ee3a5700df34321d62077e63633c575c1c954514e99da7c179d");
        assert_eq!(x448(&k, &u), r);
    }

    // RFC 7748, section 6.2: the full Diffie-Hellman exchange.
    #[test]
    fn diffie_hellman_agrees() {
        let a = SecretKey::from_bytes(h(
            "9a8f4925d1519f5775cf46b04b5800d4ee9ee8bae8bc5565d498c28dd9c9baf574a9419744897391006382a6f127ab1d9ac2d8c0a598726b",
        ));
        let b = SecretKey::from_bytes(h(
            "1c306a7ac2a0e2e0990b294470cba339e6453772b075811d8fad0d1d6927c120bb5ee8972b0d3e21374c9c921b09d1b0366f10b65173992d",
        ));

        let a_pub = a.public_key();
        let b_pub = b.public_key();
        // RFC 7748 section 6.2 expected public keys
        assert_eq!(
            a_pub.to_bytes(),
            h("9b08f7cc31b7e3e67d22d5aea121074a273bd2b83de09c63faa73d2c22c5d9bbc836647241d953d40c5b12da88120d53177f80e532c41fa0")
        );
        assert_eq!(
            b_pub.to_bytes(),
            h("3eb7a829b0cd20f5bcfc0b599b6feccf6da4627107bdb0d4f345b43027d8b972fc3e34fb4232a13ca706dcb57aec3dae07bdc1c67bf33609")
        );

        let ab = a.diffie_hellman(&b_pub);
        let ba = b.diffie_hellman(&a_pub);
        assert_eq!(ab.to_bytes(), ba.to_bytes());
        // RFC 7748 section 6.2 expected shared secret
        assert_eq!(
            ab.to_bytes(),
            h("07fff4181ac6cc95ec1c16a94a0f74d12da232ce40a77552281d282bb60c0b56fd2464c335543936521c24403085d59a449a5037514a879d")
        );
        assert!(ab.was_contributory());
    }
}
