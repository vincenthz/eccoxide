//! Divan micro-benchmarks for X25519 key agreement (RFC 7748), comparing this
//! crate (`eccoxide`) against the `cryptoxide` implementation.
//!
//! The `cryptoxide` comparison benches require the optional `cryptoxide`
//! dependency, which is pulled in by the `ed25519` feature (or can be enabled
//! directly).
//!
//! Run with:
//!
//! ```text
//! cargo bench --bench x25519 --features "x25519 ed25519"
//! cargo bench --bench x25519 --features x25519              # eccoxide only
//! ```

fn main() {
    divan::main();
}

// Two arbitrary secret scalars (RFC 7748 section 6.1).
#[cfg(feature = "x25519")]
const SK_A: [u8; 32] = [
    0x77, 0x07, 0x6d, 0x0a, 0x73, 0x18, 0xa5, 0x7d, 0x3c, 0x16, 0xc1, 0x72, 0x51, 0xb2, 0x66, 0x45,
    0xdf, 0x4c, 0x2f, 0x87, 0xeb, 0xc0, 0x99, 0x2a, 0xb1, 0x77, 0xfb, 0xa5, 0x1d, 0xb9, 0x2c, 0x2a,
];
#[cfg(feature = "x25519")]
const SK_B: [u8; 32] = [
    0x5d, 0xab, 0x08, 0x7e, 0x62, 0x4a, 0x8a, 0x4b, 0x79, 0xe1, 0x7f, 0x8b, 0x83, 0x80, 0x0e, 0xe6,
    0x6f, 0x3b, 0xb1, 0x29, 0x26, 0x18, 0xb6, 0xfd, 0x1c, 0x2f, 0x8b, 0x27, 0xff, 0x88, 0xe0, 0xeb,
];

/// Public-key derivation: `[scalar]·basepoint`.
#[cfg(feature = "x25519")]
mod derive_public {
    use super::SK_A;
    use divan::{black_box, Bencher};

    #[divan::bench]
    fn eccoxide(bencher: Bencher) {
        use ::eccoxide::protocol::x25519::x25519_base;
        bencher.bench(|| x25519_base(black_box(&SK_A)));
    }

    #[cfg(feature = "cryptoxide")]
    #[divan::bench]
    fn cryptoxide(bencher: Bencher) {
        use ::cryptoxide::curve25519::curve25519_base;
        bencher.bench(|| curve25519_base(black_box(&SK_A)));
    }
}

/// Full Diffie-Hellman: `[scalar]·peer_public`.
#[cfg(feature = "x25519")]
mod agreement {
    use super::{SK_A, SK_B};
    use divan::{black_box, Bencher};

    #[divan::bench]
    fn eccoxide(bencher: Bencher) {
        use ::eccoxide::protocol::x25519::{x25519, x25519_base};
        let peer = x25519_base(&SK_B);
        bencher.bench(|| x25519(black_box(&SK_A), black_box(&peer)));
    }

    #[cfg(feature = "cryptoxide")]
    #[divan::bench]
    fn cryptoxide(bencher: Bencher) {
        use ::cryptoxide::curve25519::{curve25519, curve25519_base};
        let peer = curve25519_base(&SK_B);
        bencher.bench(|| curve25519(black_box(&SK_A), black_box(&peer)));
    }
}
