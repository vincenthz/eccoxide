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

For now, most SEC2 curves are supported through fiat-crypto:

* p256r1, p256k1, p384r1, p521r1
* p190r1, p190k1: not particularly recommended due to size

Special cases:

* p224k1: p=5 mod 8, using alternative approach for sqrt calculation
* p224r1: p=1 mod 8, using tonelli shanks algorithm for sqrt calculation

Optionally someone can enable all SEC2 curves less than 190bits (112 to 160 bits)
using sec2-small features, but the size of those curves are too small to be used
in normal settings. Also those curves are using a generic backend using num-traits
and num-bigint, which is not particularly fast, nor secure.

Futures plans includes support of ed25519, ed448, curve9767, and other edwards curves,
and maybe other.

## FAQ

Q: Does using formally generated modules makes this crate more secure ?
A: No, while it improve basic guaranteed of correctness, it is also based on
   model that are assumed correct. It also depends on the rust/llvm compiler to
   not bring various optimisation / code change that could break some properties
   (e.g. constant time) and finally there's also lots of glue being written on top
   to provide high level usable ECC, that have been manually written.

## :warning: Disclaimer

This is not a ready-to-use in production code crate. see TODO.md.
