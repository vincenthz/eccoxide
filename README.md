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

## FAQ

Q: Does using formally generated modules makes this crate more secure ?
A: No, while it improve basic guaranteed of correctness, it is also based on
   model that are assumed correct. It also depends on the rust/llvm compiler to
   not bring various optimisation / code change that could break some properties
   (e.g. constant time) and finally there's also lots of glue being written on top
   to provide high level usable ECC, that have been manually written.

## :warning: Disclaimer

This is not a ready-to-use in production code crate. see TODO.md.
