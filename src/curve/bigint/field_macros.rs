#[doc(hidden)]
#[macro_export]
macro_rules! bigint_scalar_impl {
    ($ty: ident, $p: expr, $sz: expr, $pmod4: expr, $pp1d4: expr) => {
        #[derive(Clone)]
        pub struct $ty(num_bigint::BigUint);

        impl PartialEq for $ty {
            fn eq(&self, other: &Self) -> bool {
                &self.0 == &other.0
            }
        }

        impl std::fmt::Debug for $ty {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                let bs = self.0.to_bytes_be();
                for b in bs.iter() {
                    write!(f, "{:02x}", b)?
                }
                Ok(())
            }
        }

        impl Eq for $ty {}

        impl $ty {
            pub const SIZE_BITS: usize = $sz;
            pub const SIZE_BYTES: usize = (Self::SIZE_BITS + 7) / 8;

            /// the zero constant (additive identity)
            pub fn zero() -> Self {
                use num_traits::identities::Zero;
                Self(BigUint::zero())
            }

            /// The one constant (multiplicative identity)
            pub fn one() -> Self {
                use num_traits::identities::One;
                Self(BigUint::one())
            }

            pub fn from_u64(n: u64) -> Self {
                use num_traits::cast::FromPrimitive;
                Self(BigUint::from_u64(n).unwrap())
            }

            pub fn is_zero(&self) -> bool {
                use num_traits::identities::Zero;
                self.0.is_zero()
            }

            // there's no really negative number in Fp, but if high bit is set ...
            pub fn high_bit_set(&self) -> bool {
                //use num_traits::identities::Zero;
                use num_traits::cast::FromPrimitive;
                self.0 > ($p / BigUint::from_u64(2).unwrap())
            }

            /// Self add another Scalar
            pub fn add_assign(&mut self, other: &Self) {
                self.0 += &other.0;
                self.0 %= $p;
            }

            /// Get the multiplicative inverse
            ///
            /// Note that 0 doesn't have a multiplicative inverse
            pub fn inverse(&self) -> Option<Self> {
                /*
                use num_traits::cast::FromPrimitive;
                let pm2 = $p - BigUint::from_u64(2).unwrap();
                Some(Self(self.0.modpow(&pm2, $p)))
                */

                use num_traits::identities::Zero;
                if self.0.is_zero() {
                    None
                } else {
                    Some(Self(mod_inverse(&self.0, $p)))
                }
            }

            /// Double the field element, this is equivalent to 2*self or self+self, but can be implemented faster
            pub fn double(&self) -> Self {
                self + self
            }

            /// Compute the field element raised to a power of n, modulus p
            pub fn power(&self, n: u64) -> Self {
                Self(self.0.modpow(&n.into(), $p))
            }

            /// Compute the square root 'x' of the field element such that x*x = self
            pub fn sqrt(&self) -> Option<Self> {
                if *$pmod4 == 3 {
                    // P mod 4 == 3, then we can compute sqrt with one exponentiation with (P+1)/4
                    Some(Self(self.0.modpow(&*$pp1d4, $p)))
                } else {
                    tonelli_shanks(&self.0, $p).map(|n| Self(n))
                }
            }

            /// Initialize a new scalar from its bytes representation
            ///
            /// If the represented value overflow the field element size,
            /// then None is returned.
            pub fn from_bytes(bytes: &[u8; Self::SIZE_BYTES]) -> Option<Self> {
                let n = BigUint::from_bytes_be(bytes);
                if &n >= $p {
                    None
                } else {
                    Some(Self(n))
                }
            }

            /// Similar to from_bytes but take values from a slice.
            ///
            /// If the slice is not of the right size, then None is returned
            pub fn from_slice(slice: &[u8]) -> Option<Self> {
                if slice.len() != Self::SIZE_BYTES {
                    return None;
                }
                let n = BigUint::from_bytes_be(slice);
                if &n >= $p {
                    None
                } else {
                    Some(Self(n))
                }
            }

            /// Output the scalar bytes representation
            pub fn to_bytes(&self) -> [u8; Self::SIZE_BYTES] {
                let mut out = [0u8; Self::SIZE_BYTES];
                let bs = self.0.to_bytes_be();
                let start: usize = Self::SIZE_BYTES - bs.len();

                // skip some bytes at the beginning if necessary, act as a 0-pad
                out[start..].copy_from_slice(&bs);
                out
            }

            /// Output the scalar bytes representation to the mutable slice
            ///
            /// the slice needs to be of the correct size
            pub fn to_slice(&self, slice: &mut [u8]) {
                assert_eq!(slice.len(), Self::SIZE_BYTES);

                // TODO don't create temporary buffer
                let bytes = self.to_bytes();
                slice.copy_from_slice(&bytes[..]);
            }

            /// Initialize from a wide buffer of random data.
            ///
            /// The difference with 'from_bytes' or 'from_slice' is that it takes
            /// a random initialized buffer and used modulo operation to initialize
            /// as a field element, but due to inherent bias in modulo operation
            /// we take a double sized buffer.
            pub fn init_from_wide_bytes(random: [u8; Self::SIZE_BYTES * 2]) -> Self {
                Self(BigUint::from_bytes_be(&random) % $p)
            }
        }

        impl std::ops::Neg for $ty {
            type Output = $ty;

            fn neg(self) -> Self::Output {
                $ty($p - self.0)
            }
        }

        impl std::ops::Neg for &$ty {
            type Output = $ty;

            fn neg(self) -> Self::Output {
                $ty($p - &self.0)
            }
        }

        // ****************
        // Scalar Addition
        // ****************

        impl<'a, 'b> std::ops::Add<&'b $ty> for &'a $ty {
            type Output = $ty;

            fn add(self, other: &'b $ty) -> $ty {
                $ty((&self.0 + &other.0) % $p)
            }
        }

        impl<'a> std::ops::Add<$ty> for &'a $ty {
            type Output = $ty;

            fn add(self, other: $ty) -> $ty {
                $ty((&self.0 + &other.0) % $p)
            }
        }

        impl<'b> std::ops::Add<&'b $ty> for $ty {
            type Output = $ty;

            fn add(self, other: &'b $ty) -> $ty {
                $ty((&self.0 + &other.0) % $p)
            }
        }

        impl std::ops::Add<$ty> for $ty {
            type Output = $ty;

            fn add(self, other: $ty) -> $ty {
                $ty((&self.0 + &other.0) % $p)
            }
        }

        // *******************
        // Scalar Subtraction
        // *******************

        impl<'a, 'b> std::ops::Sub<&'b $ty> for &'a $ty {
            type Output = $ty;

            fn sub(self, other: &'b $ty) -> $ty {
                $ty((&self.0 + (-other).0) % $p)
            }
        }

        impl<'a> std::ops::Sub<$ty> for &'a $ty {
            type Output = $ty;

            fn sub(self, other: $ty) -> $ty {
                self - &other
            }
        }

        impl<'b> std::ops::Sub<&'b $ty> for $ty {
            type Output = $ty;

            fn sub(self, other: &'b $ty) -> $ty {
                &self - other
            }
        }

        impl std::ops::Sub<$ty> for $ty {
            type Output = $ty;

            fn sub(self, other: $ty) -> $ty {
                &self - &other
            }
        }

        // **********************
        // Scalar Multiplication
        // **********************

        impl<'a, 'b> std::ops::Mul<&'b $ty> for &'a $ty {
            type Output = $ty;

            fn mul(self, other: &'b $ty) -> $ty {
                $ty((&self.0 * &other.0) % $p)
            }
        }

        impl<'b> std::ops::Mul<&'b $ty> for $ty {
            type Output = $ty;

            fn mul(self, other: &'b $ty) -> $ty {
                $ty((&self.0 * &other.0) % $p)
            }
        }

        impl<'a, 'b> std::ops::Mul<$ty> for &'a $ty {
            type Output = $ty;

            fn mul(self, other: $ty) -> $ty {
                $ty((&self.0 * &other.0) % $p)
            }
        }

        impl std::ops::Mul<$ty> for $ty {
            type Output = $ty;

            fn mul(self, other: $ty) -> $ty {
                $ty((&self.0 * &other.0) % $p)
            }
        }
    };
}
