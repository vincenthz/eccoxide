//! Divan micro-benchmarks for X25519 key agreement (RFC 7748).
//!
//! Run with:
//!
//! ```text
//! cargo bench --bench x25519 --features x25519
//! ```

fn main() {
    divan::main();
}

#[cfg(feature = "x25519")]
mod x25519 {
    use divan::{black_box, Bencher};
    use eccoxide::protocol::x25519::{x25519, x25519_base, SecretKey};

    // Two arbitrary (unclamped) secret scalars.
    const SK_A: [u8; 32] = [
        0x77, 0x07, 0x6d, 0x0a, 0x73, 0x18, 0xa5, 0x7d, 0x3c, 0x16, 0xc1, 0x72, 0x51, 0xb2, 0x66,
        0x45, 0xdf, 0x4c, 0x2f, 0x87, 0xeb, 0xc0, 0x99, 0x2a, 0xb1, 0x77, 0xfb, 0xa5, 0x1d, 0xb9,
        0x2c, 0x2a,
    ];
    const SK_B: [u8; 32] = [
        0x5d, 0xab, 0x08, 0x7e, 0x62, 0x4a, 0x8a, 0x4b, 0x79, 0xe1, 0x7f, 0x8b, 0x83, 0x80, 0x0e,
        0xe6, 0x6f, 0x3b, 0xb1, 0x29, 0x26, 0x18, 0xb6, 0xfd, 0x1c, 0x2f, 0x8b, 0x27, 0xff, 0x88,
        0xe0, 0xeb,
    ];

    /// Public-key derivation: `[scalar]·basepoint`.
    #[divan::bench]
    fn derive_public(bencher: Bencher) {
        bencher.bench(|| x25519_base(black_box(&SK_A)));
    }

    /// Full Diffie-Hellman: `[scalar]·peer_public`.
    #[divan::bench]
    fn agreement(bencher: Bencher) {
        let peer = x25519_base(&SK_B);
        bencher.bench(|| x25519(black_box(&SK_A), black_box(&peer)));
    }

    /// The typed `SecretKey::diffie_hellman` wrapper, end to end.
    #[divan::bench]
    fn typed_diffie_hellman(bencher: Bencher) {
        let a = SecretKey::from_bytes(SK_A);
        let b_pub = SecretKey::from_bytes(SK_B).public_key();
        bencher.bench(|| black_box(&a).diffie_hellman(black_box(&b_pub)));
    }
}
