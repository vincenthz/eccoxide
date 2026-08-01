# Changelog

## Unreleased

## 0.4.3 - 2026-08-01

- **SM2**: the recommended 256-bit prime-field curve of the Chinese SM2 standard
  (GB/T 32918 / GM/T 0003), a short Weierstrass curve with `a = -3` wired like
  the SEC2 `r1` curves, including a fixed-base comb table.
- **The BLS12-381 pairing is roughly 24x faster** (~28.7 ms to ~1.2 ms):
- `PointAffine::decompress` is now constant time.
- Updated `cryptoxide` to 0.6.

## 0.4.2 - 2026-07-18

- **BLS12-381**: the base field `Fp`, the scalar field, the full `Fp2`/`Fp6`/`Fp12`
  extension tower, the G1 and G2 groups, and the optimal-ate pairing
  `e: G1 x G2 -> Fp12`, plus the fiat-crypto extraction for both fields
  (64-bit only).
- **Jubjub**: the twisted Edwards curve defined over the BLS12-381 scalar field.

## 0.4.1 - 2026-07-11

- **ECDSA** signatures over the SEC2 Weierstrass curves.
- **ristretto255**, built on top of curve25519.

## 0.4.0 - 2026-06-22

- **curve25519 / edwards25519**, with the **X25519** and **Ed25519** protocols,
  the Edwards and Montgomery form abstractions, and a generator table.
- **curve448** and the **X448** protocol.
- Fixed-base (generator) comb precomputation tables and `mul_base`, enabled by
  default through the new `table` feature.
- Variable-time scalar multiplication `mul_vartime` using wNAF (window non-adjacent form)
- Bernstein-Yang safegcd modular inversion.
- `init_from_wide_bytes` support.
- Ed25519 keypair signing
- Benchmarks for the SEC2 curves, curve25519 and x448.
- Point addition on several curves now uses the complete addition formulas
  instead of branching on equality to choose between add and double.

## 0.3.1 - 2021-10-24

- Support for the fiat-crypto `unsaturated_solinas` representation.

## 0.3.0 - 2020-10-29

- **secp192r1**, **secp192k1**, **secp224r1** and **secp224k1**

## 0.2.0 - 2020-10-01

- **secp256r1**, **secp384r1** and **secp521r1**
- Point negation, exposed curve parameters, and a `Sign` type.
- Sage-generated known-answer tests.
- Documentation.

## 0.1.0

- Initial version: **secp256k1**
