## General

features:

* [ ] NonZeroFieldElement / NonZeroScalar types
* [ ] constant-time Result/Either (CtResult, CtEither)
* [ ] assign-style operators (add_assign / sub_assign / mul_assign)
* [ ] Scalar -> FieldElement conversion
* [ ] generic hash-to-curve (point)
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
* [ ] hash-to-curve (RFC 9380) for G1 and G2: SSWU map-to-curve + isogeny maps,
  expand_message_xmd, and the encode/hash variants
* [ ] subgroup membership checks (is-in-G1 / is-in-G2) and cofactor clearing
* [ ] point serialization in the standard zcash/IETF compressed & uncompressed
  formats (the infinity / compression / sort flag bits)
* [ ] BLS signatures (min-pk and min-sig variants), signature aggregation and
  proof-of-possession
* [ ] multi-pairing / pairing product: accumulate several Miller loops then a
  single final exponentiation

optimisation (pairing is currently correctness-first, not fast):
* [ ] projective coordinates + sparse line functions (mul_by_014) in the Miller loop
* [ ] Frobenius-based final exponentiation (easy + hard part) instead of the single
  square-and-multiply by (p^12-1)/r ; precompute the Frobenius coefficients
* [ ] cyclotomic squaring in the final exponentiation
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
* [ ] cached ("niels") point representation (y-x, y+x, 2d*T) so the hot addition does
  not recompute those each call
* [ ] mixed-coordinate addition (assume Z2 = 1) for the affine comb-table points in
  `mul_base` and in the table build, saving field multiplications per step
* [ ] larger or signed-digit comb windows for `mul_base` to cut additions further
* [ ] X25519: the Montgomery ladder is already standard x-only add-and-double, so the
  remaining lever is the field multiply — benchmark the saturated 4-limb
  (dettman) fiat backend against the current 5-limb unsaturated-solinas one
* [ ] the same variable-base speedups apply to every other Edwards/Weierstrass
  `Point::scale` (jubjub, edwards448, sec2)

## ed25519 (performance)

* [ ] verify currently does a constant-time `mul_base(S)` plus a full constant-time
  variable-base `scale(A)` and then compares; verification is over public data,
  so replace it with simultaneous double-scalar multiplication for
  [S]B - [k]A using wNAF and a precomputed table for B (variable-time)
* [ ] add a variable-time (w)NAF variable-base multiplication as the shared building
  block for the above
* [ ] batch verification: verify n signatures with one multi-scalar multiplication
  over a random linear combination of the equations
* [ ] cache a niels-form precomputation of a repeatedly-used public key A across
  verifications
* conformance / strict validation criteria ("it's 255:19 AM" blog bost - ed25519-zebra)
