//! Divan micro-benchmarks for Ed25519 signatures (RFC 8032), comparing this
//! crate (`eccoxide`) against the `cryptoxide` implementation.
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
const SEED: [u8; 32] = [
    0x9d, 0x61, 0xb1, 0x9d, 0xef, 0xfd, 0x5a, 0x60, 0xba, 0x84, 0x4a, 0xf4, 0x92, 0xec, 0x2c, 0xc4,
    0x44, 0x49, 0xc5, 0x69, 0x7b, 0x32, 0x69, 0x19, 0x70, 0x3b, 0xac, 0x03, 0x1c, 0xae, 0x7f, 0x60,
];

#[cfg(feature = "ed25519")]
const MESSAGE: &[u8] = b"the quick brown fox jumps over the lazy dog";

/// Key generation: SHA-512 of the seed plus one fixed-base scalar mult.
#[cfg(feature = "ed25519")]
mod keygen {
    use super::SEED;
    use divan::{black_box, Bencher};

    #[divan::bench]
    fn eccoxide(bencher: Bencher) {
        use ::eccoxide::protocol::ed25519::Keypair;
        bencher.bench(|| Keypair::from_seed(black_box(SEED)));
    }

    #[divan::bench]
    fn cryptoxide(bencher: Bencher) {
        bencher.bench(|| ::cryptoxide::ed25519::keypair(black_box(&SEED)));
    }
}

/// Signing a message.
#[cfg(feature = "ed25519")]
mod sign {
    use super::{MESSAGE, SEED};
    use divan::{black_box, Bencher};

    #[divan::bench]
    fn eccoxide(bencher: Bencher) {
        use ::eccoxide::protocol::ed25519::SecretKey;
        let sk = SecretKey::from_bytes(SEED);
        bencher.bench(|| black_box(&sk).sign(black_box(MESSAGE)));
    }

    #[divan::bench]
    fn cryptoxide(bencher: Bencher) {
        let (keypair, _public) = ::cryptoxide::ed25519::keypair(&SEED);
        bencher.bench(|| ::cryptoxide::ed25519::signature(black_box(MESSAGE), black_box(&keypair)));
    }
}

/// Verifying a signature.
#[cfg(feature = "ed25519")]
mod verify {
    use super::{MESSAGE, SEED};
    use divan::{black_box, Bencher};

    #[divan::bench]
    fn eccoxide(bencher: Bencher) {
        use ::eccoxide::protocol::ed25519::SecretKey;
        let sk = SecretKey::from_bytes(SEED);
        let pk = sk.public_key();
        let sig = sk.sign(MESSAGE);
        bencher.bench(|| black_box(&pk).verify(black_box(MESSAGE), black_box(&sig)));
    }

    #[divan::bench]
    fn cryptoxide(bencher: Bencher) {
        let (keypair, public) = ::cryptoxide::ed25519::keypair(&SEED);
        let sig = ::cryptoxide::ed25519::signature(MESSAGE, &keypair);
        bencher
            .bench(|| ::cryptoxide::ed25519::verify(black_box(MESSAGE), black_box(&public), black_box(&sig)));
    }
}
