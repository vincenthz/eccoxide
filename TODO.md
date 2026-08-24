## General

features:

* [ ] NonZeroFieldElement / NonZeroScalar types
* [ ] constant-time Result/Either (CtResult, CtEither)
* [ ] assign-style operators (add_assign / sub_assign / mul_assign)
* [ ] Scalar -> FieldElement conversion
* [ ] generic hash-to-curve (point), factoring out the BLS12-381 RFC 9380
  implementation (SSWU + isogeny + expand_message) for the other curves
* [ ] RFC 6979 deterministic nonce generation for ECDSA
* [ ] Schnorr signatures
* [ ] Ed448 signatures on curve448
* [ ] Pairing related protocols

optimisation:
* [ ] fill in the remaining fiat backends: generated sqrt & inverse addition-chains,
  a macro to write addition chains, and constant-time variants

tooling:
* [ ] audit functions for constant-time-ness

## BLS12-381

features:
* [ ] Generalize Hash algorithms to be able to pass arbitrary hash functions that fullfill
  the output size requirement.
* [ ] `expand_message_xof` (SHAKE128) alongside the `_XMD:SHA-256_` suites, for
  applications that ask for the XOF variants of RFC 9380
* [ ] GT membership check (an element of the cyclotomic subgroup with `f^p == f^x`),
  the counterpart of the G1/G2 `is_in_subgroup` checks for pairing outputs
* [ ] BLS signatures (min-pk and min-sig variants), signature aggregation and
  proof-of-possession

optimisation
* [ ] drop the inversion `map_to_curve` ends with (`x = x / tv4`, step 25 of
  RFC 9380 appendix F.2) by handing `(x_num, x_den)` to the isogeny, which is
  already evaluated projectively: a homogeneous `poly_eval` would need the
  denominator powers but no `Fp`/`Fp2` inverse. Two inversions per
  `hash_to_curve`, ~24us each, out of ~156us for G1 and ~625us for G2
* [ ] windowed / NAF exponentiation for the 126-bit lambda3 step, whose 48 set bits
  cost ~47 of the hard part's ~68 Fp12 multiplications
* [ ] compressed cyclotomic squaring (Karabina) for the exp-by-x chains
* [ ] GLV endomorphism on G1 (and psi on G2) for faster scalar multiplication
* [ ] multi-scalar multiplication (Pippenger / wNAF)
* [ ] hand-written Fermat addition chains for Fp / Fr / Fp2 inverse (currently
  safegcd), and an Fp2 sqrt/inverse without falling back to the base field

## Jubjub

features:
* [ ] Montgomery-form birational map (Montgomery A = 40962) and x-only Montgomery
  ladder, mirroring curve25519
* [ ] packed 32-byte compressed point encoding (little-endian y with the x sign bit),
  as used by zcash
* [ ] hash-to-curve (Elligator 2 via the Montgomery form, or the twisted-Edwards map)
* [ ] full-order (cofactor-8) point handling: cofactor multiplication and prime-order
  subgroup membership checks

optimisation:
* [ ] constant-time square root in bls12_381::Scalar (Fr) — the current
  Tonelli-Shanks used by point decompression is variable-time
* [ ] Scalar inverse via a hand-written addition chain instead of safegcd

refactor:
* [ ] factor the duplicated edwards25519 / jubjub `Point` code into a generic
  complete twisted-Edwards `Point<C: EdwardsCurveAM1>` shared by both

## curve25519 / edwards25519 (performance)

closing the gap vs the cryptoxide comparison benches (benches/curve25519, x25519):

* [ ] variable-base scalar multiplication (`Point::scale` / `scale_bytes`) is a naive
  constant-time bit-by-bit double-and-add (~256 doublings + 256 masked additions
  per scalar); switch to a fixed 4-bit window with a constant-time point table
  (~256 D + 64 A), or a signed-window / (w)NAF method for the variable-time paths
* [ ] the width-`WNAF_BASE_W` generator table of
  `double_scalar_mul_base_vartime` is still `CachedPoint`, so every base
  addition pays a multiplication by a `z` that is one: move it to
  `CachedPointAffine` / `add_cached_affine` as `mul_base` now does
* [ ] `add_wnaf_multiple` builds `entry.negate()` for a negative digit, four
  field-element clones per addition; negate inside the formula instead
* [ ] X25519: the Montgomery ladder is already standard x-only add-and-double, so the
  remaining lever is the field multiply — benchmark the saturated 4-limb
  (dettman) fiat backend against the current 5-limb unsaturated-solinas one
* [ ] the same variable-base speedups apply to every other Edwards/Weierstrass
  `Point::scale` (jubjub, edwards448, sec2)

## ed25519 (performance)

* [ ] batch verification: verify n signatures with one multi-scalar multiplication
  over a random linear combination of the equations
* [ ] cache a niels-form precomputation of a repeatedly-used public key A across
  verifications
* conformance / strict validation criteria ("it's 255:19 AM" blog bost - ed25519-zebra)
