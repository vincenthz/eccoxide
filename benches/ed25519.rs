//! Divan micro-benchmarks for Ed25519 signatures (RFC 8032).
//!
//! Run with:
//!
//! ```text
//! cargo bench --bench ed25519 --features ed25519
//! ```

fn main() {
    divan::main();
}

#[cfg(feature = "ed25519")]
mod ed25519 {
    use divan::{black_box, Bencher};
    use eccoxide::protocol::ed25519::{Keypair, PublicKey, SecretKey, Signature};

    const SEED: [u8; 32] = [
        0x9d, 0x61, 0xb1, 0x9d, 0xef, 0xfd, 0x5a, 0x60, 0xba, 0x84, 0x4a, 0xf4, 0x92, 0xec, 0x2c,
        0xc4, 0x44, 0x49, 0xc5, 0x69, 0x7b, 0x32, 0x69, 0x19, 0x70, 0x3b, 0xac, 0x03, 0x1c, 0xae,
        0x7f, 0x60,
    ];

    const MESSAGE: &[u8] = b"the quick brown fox jumps over the lazy dog";

    /// Key generation: SHA-512 of the seed plus one fixed-base scalar mult.
    #[divan::bench]
    fn keygen(bencher: Bencher) {
        bencher.bench(|| Keypair::from_seed(black_box(SEED)));
    }

    /// Public-key derivation from an existing secret key.
    #[divan::bench]
    fn public_key(bencher: Bencher) {
        let sk = SecretKey::from_bytes(SEED);
        bencher.bench(|| black_box(&sk).public_key());
    }

    /// Signing a message.
    #[divan::bench]
    fn sign(bencher: Bencher) {
        let sk = SecretKey::from_bytes(SEED);
        bencher.bench(|| black_box(&sk).sign(black_box(MESSAGE)));
    }

    /// Verifying a signature.
    #[divan::bench]
    fn verify(bencher: Bencher) {
        let sk = SecretKey::from_bytes(SEED);
        let pk: PublicKey = sk.public_key();
        let sig: Signature = sk.sign(MESSAGE);
        bencher.bench(|| black_box(&pk).verify(black_box(MESSAGE), black_box(&sig)));
    }
}
