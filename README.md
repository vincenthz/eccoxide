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

Most SEC2 curves are supported through fiat-crypto:

* p256r1, p256k1, p384r1, p521r1
* p192r1, p192k1: not particularly recommended due to size
* p224k1: p=5 mod 8, using alternative approach for sqrt calculation
* p224r1: p=1 mod 8, using (non constant time) tonelli shanks algorithm for sqrt calculation

Optionally someone can enable all SEC2 curves less than 192bits (112 to 160 bits)
using the `sec2-small` feature, but the size of those curves are too small to be used
in normal settings. Also those curves are using a generic backend using num-traits
and num-bigint, which is not particularly fast, nor secure.

### Edwards / Montgomery

* curve25519: Montgomery Curve25519 and Twisted Edwards Ed25519
* curve448: The Goldilocks Curve448 (Montgomery) and edwards448 (Edwards) curves

## Protocols

Higher-level protocols are built on top of the curves above (in the `protocol` module),
each behind its own cargo feature:

* `x25519`: X25519 Diffie-Hellman key agreement (RFC 7748), on Curve25519
* `ed25519`: Ed25519 digital signatures (RFC 8032), on edwards25519
* `x448`: X448 Diffie-Hellman key agreement (RFC 7748), on Curve448

## Features

* `sec2` (default): all SEC2 curves from 192 bits and up
* `sec2-small`: the smaller SEC2 curves (112 to 160 bits), via the generic bigint backend
* `table` (default): embed fixed-base precomputation tables so `Point::mul_base` uses a
  constant-time comb (~4x faster); adds static data to the binary
* `curve25519` (default), `curve448`: the Edwards/Montgomery curves
* `x25519`, `ed25519`, `x448`: the protocols above
* individual SEC2 curves (e.g. `p256r1`) can be enabled one at a time

## Future plans

Future plans include Ed448 signatures, curve9767, hash-to-curve, ECDSA on the
Weierstrass curves, and other curves.

## FAQ

Q: Does using formally generated modules makes this crate more secure ?
A: No, while it improve basic guaranteed of correctness, it is also based on
   model that are assumed correct. It also depends on the rust/llvm compiler to
   not bring various optimisation / code change that could break some properties
   (e.g. constant time) and finally there's also lots of glue being written on top
   to provide high level usable ECC, that have been manually written.

## :warning: Disclaimer

This is not a ready-to-use in production code crate. see TODO.md.
