//! Divan micro-benchmarks for Ed25519 signatures (RFC 8032), comparing this
//! crate (`eccoxide`) against the `cryptoxide` implementation, and the plain
//! verification against the precomputed one (`PrecomputedPublicKey`).
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

    // Both sides sign from a precomputed keypair (cached public key), so this is
    // an apples-to-apples comparison of the per-signature cost.
    #[divan::bench]
    fn eccoxide(bencher: Bencher) {
        use ::eccoxide::protocol::ed25519::Keypair;
        let kp = Keypair::from_seed(SEED);
        bencher.bench(|| black_box(&kp).sign(black_box(MESSAGE)));
    }

    #[divan::bench]
    fn cryptoxide(bencher: Bencher) {
        let (keypair, _public) = ::cryptoxide::ed25519::keypair(&SEED);
        bencher.bench(|| ::cryptoxide::ed25519::signature(black_box(MESSAGE), black_box(&keypair)));
    }
}

/// Verifying a signature.
///
/// `eccoxide_precomputed` verifies with the per-key precomputation already
/// done, which is what verifying many signatures under one key amortises down
/// to; `precompute` below is what it costs to get there.
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
    fn eccoxide_precomputed(bencher: Bencher) {
        use ::eccoxide::protocol::ed25519::SecretKey;
        let sk = SecretKey::from_bytes(SEED);
        let pk = sk.public_key().precompute().expect("a key we just derived");
        let sig = sk.sign(MESSAGE);
        bencher.bench(|| black_box(&pk).verify(black_box(MESSAGE), black_box(&sig)));
    }

    #[divan::bench]
    fn cryptoxide(bencher: Bencher) {
        let (keypair, public) = ::cryptoxide::ed25519::keypair(&SEED);
        let sig = ::cryptoxide::ed25519::signature(MESSAGE, &keypair);
        bencher.bench(|| {
            ::cryptoxide::ed25519::verify(black_box(MESSAGE), black_box(&public), black_box(&sig))
        });
    }
}

/// Building the precomputation: one point decompression (a field square root)
/// and one table of odd multiples (a single field inversion, shared).
#[cfg(feature = "ed25519")]
mod precompute {
    use super::SEED;
    use divan::{black_box, Bencher};

    #[divan::bench]
    fn eccoxide(bencher: Bencher) {
        use ::eccoxide::protocol::ed25519::SecretKey;
        let pk = SecretKey::from_bytes(SEED).public_key();
        bencher.bench(|| black_box(&pk).precompute().unwrap());
    }
}

/// Verifying `n` signatures under a single key, the precomputation counted in,
/// which is where the two paths can be compared honestly: `n = 1` is what the
/// precomputation costs on top of a single verification, and the point where
/// the curves cross is how many signatures a key has to verify to be worth
/// precomputing.
#[cfg(feature = "ed25519")]
mod verify_n {
    use super::{MESSAGE, SEED};
    use divan::{black_box, Bencher};
    use eccoxide::protocol::ed25519::{SecretKey, Signature};

    const COUNTS: &[usize] = &[1, 2, 4, 16, 64];

    /// `n` signatures over `n` distinct messages, under one key.
    fn signatures(n: usize) -> (SecretKey, Vec<Vec<u8>>, Vec<Signature>) {
        let sk = SecretKey::from_bytes(SEED);
        let messages: Vec<Vec<u8>> = (0..n)
            .map(|i| {
                let mut m = MESSAGE.to_vec();
                m.push(i as u8);
                m
            })
            .collect();
        let signatures = messages.iter().map(|m| sk.sign(m)).collect();
        (sk, messages, signatures)
    }

    #[divan::bench(args = COUNTS)]
    fn eccoxide(bencher: Bencher, n: usize) {
        let (sk, messages, signatures) = signatures(n);
        let pk = sk.public_key();
        bencher.bench(|| {
            let pk = black_box(&pk);
            let mut ok = true;
            for (m, sig) in messages.iter().zip(signatures.iter()) {
                ok &= pk.verify(black_box(m), black_box(sig));
            }
            ok
        });
    }

    #[divan::bench(args = COUNTS)]
    fn eccoxide_precomputed(bencher: Bencher, n: usize) {
        let (sk, messages, signatures) = signatures(n);
        let pk = sk.public_key();
        bencher.bench(|| {
            // the precomputation is inside the measured loop: this is the whole
            // cost of verifying `n` signatures starting from the encoded key
            let pk = black_box(&pk).precompute().unwrap();
            let mut ok = true;
            for (m, sig) in messages.iter().zip(signatures.iter()) {
                ok &= pk.verify(black_box(m), black_box(sig));
            }
            ok
        });
    }

    #[divan::bench(args = COUNTS)]
    fn cryptoxide(bencher: Bencher, n: usize) {
        let (_, messages, _) = signatures(n);
        let (keypair, public) = ::cryptoxide::ed25519::keypair(&SEED);
        let sigs: Vec<[u8; 64]> = messages
            .iter()
            .map(|m| ::cryptoxide::ed25519::signature(m, &keypair))
            .collect();
        bencher.bench(|| {
            let mut ok = true;
            for (m, sig) in messages.iter().zip(sigs.iter()) {
                ok &=
                    ::cryptoxide::ed25519::verify(black_box(m), black_box(&public), black_box(sig));
            }
            ok
        });
    }
}
