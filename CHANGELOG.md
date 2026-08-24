# Changelog

## Unreleased

- **curve25519**: `Point::mul_base` is faster. The
  fixed-base comb reads *signed* 4-bit digits, so a window holds the eight
  magnitudes `1..=8` instead of fifteen digits and a constant-time lookup walks
  half as much table; the entries are stored in the "niels" form the addition
  formula reads (`y - x`, `y + x`, `2d*x*y`) rather than as affine `(x, y)`, a
  quarter less again and nothing to recompute per addition; and since their `z`
  is one, the addition is the mixed-coordinate one, seven field
  multiplications instead of nine. The embedded table is 12KB smaller for it.
- **All curves**: `init_from_wide_bytes_be` / `_le`, which reduce a
  double-width buffer into a field or scalar, evaluate the wide integer in
  half-field-wide digits instead of one byte at a time.
- **BLS12-381**: optimise `Fp6` squaring.
- **curve25519**: Ed25519 verification is faster. two scalar
  multiplications are now a single variable-time double-scalar multiplication,
  `Point::double_scalar_mul_base_vartime`
- **curve25519**: `Point::decompress` is about twice as fast: it takes the
  square root of `(y² - 1) / (d·y² + 1)` as a quotient, one exponentiation over
  the shared `u·v³·(u·v⁷)^((p-5)/8)` candidate, instead of inverting the
  denominator and then rooting, which cost two. The RFC 9496 `SQRT_RATIO_M1` of
  ristretto255 now shares that candidate.
- **Ed25519**: verification is faster: Decoding a point no longer recovers `x`
  to reject the non-canonical `(x = 0, sign = 1)` encodings: on the curve `x` is
  zero exactly for `y = ±1`, so the test is on `y` and the field inversion it
  used to cost is gone.
- **curve25519**: `Point::scale_vartime` is the width-5 wNAF variable-time
  scalar multiplication on its own, about twice as fast as the constant-time
  `Point::scale`, and is what `CurveGroup::mul_vartime` now uses.

## 0.5.0 - 2026-08-18

- **BLS12-381**: `is_in_subgroup` on the `Point` / `PointAffine` types of `g1`
  and `g2` tells the prime-order-`r` subgroup apart from the curve it sits in,
  which the curve equation alone does not. Both are the fast endomorphism checks
  of <https://eprint.iacr.org/2021/1130>: `σ(P) = [-x²]P` for G1 (`σ` being the
  cube-root-of-unity endomorphism) and `ψ(Q) = [x]Q` for G2 (`ψ` being
  untwist-Frobenius-twist), which cost a couple of multiplications by the 64-bit
  seed instead of one by `r`.
- **BLS12-381**: hashing to the groups, behind the new
  `bls12-381-hash-to-curve` feature (on by default): `Point::hash_to_curve` and
  `Point::encode_to_curve` on `g1` and `g2` implement the four RFC 9380 suites
  `BLS12381G{1,2}_XMD:SHA-256_SSWU_{RO,NU}_`, and the `hash_to_curve` module
  exposes `expand_message_xmd` on its own. The RFC's test vectors are checked
  against, down to the intermediate field elements and pre-cofactor points.
- **BLS12-381**: `clear_cofactor` on the `Point` / `PointAffine` types of `g1`
  and `g2` sends an arbitrary point of the curve (resp. the twist) into the
  prime-order subgroup, the step hash-to-curve ends with. G1 multiplies by
  `1 - x`, G2 uses the endomorphism chain of
  <https://eprint.iacr.org/2017/419>; both are the `h_eff` of RFC 9380 section
  8.8, so neither pays for a multiplication by the cofactor itself (126 and 507
  bits).
- **BLS12-381**: `from_compressed` / `from_uncompressed` now reject points
  outside the prime-order subgroup; `from_compressed_oncurve_only` /
  `from_uncompressed_oncurve_only` keep decoding any point of the curve, for when
  membership is known or checked elsewhere.
- **BLS12-381**: the two halves of the pairing are now public API: `miller_loop`
  returns a `MillerLoopResult` (which multiplies) and
  `MillerLoopResult::final_exponentiation` completes it.
- **BLS12-381**: G1 and G2 points serialize to the standard zcash/IETF
  compressed and uncompressed encodings, through `to_compressed` /
  `from_compressed` / `to_uncompressed` / `from_uncompressed` on `Point` and
  `PointAffine`. Encoding is constant time, the identity included. Decoding
  validates the flag bits, the canonicity of each coordinate and the curve
  equation;
- **Breaking**: `Fp2::to_bytes` and `Fp2::from_bytes_unchecked` now use the
  `c1 || c0` component order, the imaginary part first, instead of `c0 || c1`.
  This is the order the standard BLS12-381 encodings use, so `Fp2` bytes are now
  directly comparable with published constants and test vectors. The change is
  silent for code that only round-trips through the crate, but bytes persisted
  by an earlier version decode to the conjugate-swapped element and have to be
  swapped on read. The BLS12-381 parameters and the G2 comb table were re-encoded
  accordingly; no other API is affected.
- `Fp2` gained `from_bytes` (the canonicity-checking counterpart of
  `from_bytes_unchecked`), plus the `from_slice` / `to_slice` pair and the
  constant-time `from_bytes_ct`, mirroring what `Fp` already provides.
- Every fiat field element gained `from_bytes_ct`, the constant-time
  counterpart of `from_bytes` in the type's default byte order, and affine
  points gained `from_coordinate_ct`, which checks the curve equation without
  branching on the coordinates.

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
