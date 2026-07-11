(this is a not exhaustive list)

* NonZeroFieldElement, NonZeroScalar
* Constant time Result/Either (CtResult, CtEither)
* add assign{add,sub,mul}
* add all other fiat implementation
  * generated sqrt & inverse "addition-chain"
  * macro to write addition chain
  * constantness
* scaling functions
* Scalar to FieldElement
* NonZeroScalar to NonZeroFieldElement
* audit function for CT
* "hash"-to-curve (point)
* RFC 6979 deterministic nonce generation for ECDSA
* Ed448 signatures on curve448
* special weirstrass curves : A=0, A=-3
