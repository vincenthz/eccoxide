# ECCoxide -- Rust Elliptic Curve Cryptography

General elliptic curve cryptography

## Design space

The aim of this crate is to provide all necessary elliptic curve cryptographic to experiments
and write new protocol, based on reasonably guaranteed based modules.

The primitives needed for basic arithmetic in finite field is provided
by formally written/generated rust modules [fiat-crypto](https://github.com/mit-plv/fiat-crypto),
which aim to provide correct, secure and constant time functions to implement those finite field.

Some other parts have been generated algorithmically, so as increase the number
of supported features and reduce the work needed to provide complete set of
features for wide-variety of curves, but with no guaranteed of being the fastest

The package rely on the following priorities list:

> make it work, then make it secure, then make it fast

Finally we rely on other arithmetic tools (e.g. sage and magma) to provides some further guarantees
on the values expected.

## Curves

### Short Weierstrass (SEC2)

Most SEC2 curves are supported:

* p256r1, p256k1, p384r1, p521r1
* p192r1, p192k1: not particularly recommended due to size
* p224k1: p=5 mod 8, using alternative approach for sqrt calculation
* p224r1: p=1 mod 8, using (non constant time) tonelli shanks algorithm for sqrt calculation

The `sec2-small` feature enables p112r1, p128r1, p160k1, p160r1, and p160r2.
These curves are too small for normal cryptographic use.
The 160-bit curves use 21-byte (161-bit) scalars and 20-byte field elements.

Note p112r2 and p128r2 are not supported due to non prime order curves (cofactor 4)

Also available:

* Pairing friendly BLS12-381
* SM2: Chinese standard GB/T 32918 / GM.T 0003

### Edwards / Montgomery

* curve25519: Montgomery Curve25519 and Twisted Edwards Ed25519
* curve448: The Goldilocks Curve448 (Montgomery) and edwards448 (Edwards) curves
* jubjub: the twisted Edwards curve `-x^2 + y^2 = 1 + d*x^2*y^2` over the
  BLS12-381 scalar field (the `bls12-381` feature provides the base field).

## Protocols

Higher-level protocols are built on top of the curves above (in the `protocol` module),
each behind its own cargo feature:

* `x25519`: X25519 Diffie-Hellman key agreement (RFC 7748), on Curve25519
* `ed25519`: Ed25519 digital signatures (RFC 8032), on edwards25519
* `x448`: X448 Diffie-Hellman key agreement (RFC 7748), on Curve448
* `ristretto255`: the ristretto255 prime-order group (RFC 9496), on edwards25519,
  with canonical encoding/decoding and a uniform-bytes one-way map
* `ecdsa`: ECDSA digital signatures (SEC1), on any of the SEC2 weierstrass curves;
  sign/verify are generic over an `EcdsaOperations` trait bundling the curve
  group (a `CurveGroup`) with the message-to-scalar hashing, with SHA-2
  instantiations provided for every curve and any other hash pluggable

## Features

Curves:

* `p256k1` (default): P256K1 curve, also known as K256, widely known through bitcoin / ethereum.
* `p256r1` (default): P256R1 curve, widely used in HSMs.
* `p384r1` (default): P384R1 curve, often found in TLS.
* `p521r1` (default): P512R1 curve.
* `p224k1`: P224K1 curve
* `p224r1`: P256R1 curve, non recommended for use, non constant SQRT due to characteristic of prime field (p=1 mod 8).
* `p192k1`, `p192r1`: 192bits curves, not recommended due to size, but sometimes still used
* `p112r1`, `p128r1`, `p160k1`, `p160r1`, `p160r2`: Smaller SEC2 curves, not recommended to their small size that are under typical security margin
* `curve25519` (default), `curve448`: the Edwards/Montgomery curves
* `bls12-381` (default): the BLS12-381 pairing-friendly curve — G1, G2, the
  `Fp2`/`Fp6`/`Fp12` tower and the optimal-ate pairing
* `bls12-381-hash-to-curve` (default): hashing to G1/G2 per RFC 9380
  (pulls `cryptoxide` for SHA-256)
* `jubjub`: the Jubjub twisted Edwards curve over BLS12-381 `Fr`
* `ristretto255`

Or as group of curves:

* `sec2` (default): all SEC2 curves from 224 bits and up
* `sec2-small`: the smaller SEC2 curves (112 to 192 bits), generally not a good idea to enable

Optimisation features:

* `table` (default): embed fixed-base precomputation tables so `Point::mul_base` uses a
  constant-time comb (~4x faster); adds static data to the binary

Protocols:

* `protocols`: to enable all protocols in this crates
* `x25519`: enable X25519 (DH)
* `ed25519`: enable Ed25519 signature
* `x448`: enable X448 (DH)
* `ecdsa`: enable ECDSA for many SEC2 curves

## Future plans

Future plans include Ed448 signatures, curve9767, hash-to-curve, and other
curves.

## FAQ

Q: Does using formally generated modules makes this crate more secure ?
A: No, while it improve basic guaranteed of correctness, it is also based on
   model that are assumed correct. It also depends on the rust/llvm compiler to
   not bring various optimisation / code change that could break some properties
   (e.g. constant time) and finally there's also lots of glue being written on top
   to provide high level usable ECC, that have been manually written.

## :warning: Disclaimer

This is not a ready-to-use in production code crate. see TODO.md.
