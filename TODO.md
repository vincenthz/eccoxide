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
* init from wide binary : remove bias
  * need barrett reduction 
* "hash"-to-curve (point)
* add ECDH/ECDSA
* fence bigint implementation behind a rust package flag
* special weirstrass curves : A=0, A=-3
* non weirstrass curves
