#[doc(hidden)]
#[macro_export]
macro_rules! scalar_impl {
    ($p: expr) => {
        #[derive(Debug, Clone)]
        pub struct Scalar(num_bigint::BigUint);

        impl PartialEq for Scalar {
            fn eq(&self, other: &Self) -> bool {
                &self.0 == &other.0
            }
        }

        impl Eq for Scalar {}

        impl Scalar {
            /// the zero constant (additive identity)
            pub fn zero() -> Self {
                use num_traits::identities::Zero;
                Scalar(BigUint::zero())
            }

            /// The one constant (multiplicative identity)
            pub fn one() -> Self {
                use num_traits::identities::One;
                Scalar(BigUint::one())
            }

            /// Self add another Scalar
            pub fn add_assign(&mut self, other: &Self) {
                self.0 += &other.0;
                self.0 %= $p;
            }

            /// Negate the scalar
            pub fn negate(&self) -> Self {
                Scalar(($p - &self.0) % $p)
            }

            /// Get the multiplicative inverse
            ///
            /// Note that 0 doesn't have a multiplicative inverse
            pub fn inverse(&self) -> Option<Self> {
                use num_traits::identities::Zero;
                if self.0.is_zero() {
                    None
                } else {
                    Some(Scalar(mod_inverse(&self.0, $p)))
                }
            }

            pub fn double(&self) -> Self {
                self * self
            }

            pub fn power(&self, n: u64) -> Self {
                Scalar(self.0.modpow(&n.into(), $p))
            }

            pub fn from_bytes(slice: &[u8]) -> Option<Self> {
                let n = BigUint::from_bytes_be(&slice);
                if &n >= $p {
                    None
                } else {
                    Some(Scalar(n))
                }
            }

            pub fn to_bytes(&self) -> [u8; $sz] {
                let mut out = [0u8; $sz];
                let bytes: usize = ((self.0.bits() + 7) >> 3) as usize;
                let start: usize = $sz - bytes;

                let bs = self.0.to_bytes_be();
                out[start..].copy_from_slice(&bs);
                out
            }
        }

        // ****************
        // Scalar Addition
        // ****************

        impl<'a, 'b> std::ops::Add<&'b Scalar> for &'a Scalar {
            type Output = Scalar;

            fn add(self, other: &'b Scalar) -> Scalar {
                Scalar((&self.0 + &other.0) % $p)
            }
        }

        impl<'a> std::ops::Add<Scalar> for &'a Scalar {
            type Output = Scalar;

            fn add(self, other: Scalar) -> Scalar {
                Scalar((&self.0 + &other.0) % $p)
            }
        }

        impl<'b> std::ops::Add<&'b Scalar> for Scalar {
            type Output = Scalar;

            fn add(self, other: &'b Scalar) -> Scalar {
                Scalar((&self.0 + &other.0) % $p)
            }
        }

        impl std::ops::Add<Scalar> for Scalar {
            type Output = Scalar;

            fn add(self, other: Scalar) -> Scalar {
                Scalar((&self.0 + &other.0) % $p)
            }
        }

        // *******************
        // Scalar Subtraction
        // *******************

        impl<'a, 'b> std::ops::Sub<&'b Scalar> for &'a Scalar {
            type Output = Scalar;

            fn sub(self, other: &'b Scalar) -> Scalar {
                Scalar((&self.0 + &other.negate().0) % $p)
            }
        }

        impl<'a> std::ops::Sub<Scalar> for &'a Scalar {
            type Output = Scalar;

            fn sub(self, other: Scalar) -> Scalar {
                self - &other
            }
        }

        impl<'b> std::ops::Sub<&'b Scalar> for Scalar {
            type Output = Scalar;

            fn sub(self, other: &'b Scalar) -> Scalar {
                &self - other
            }
        }

        impl std::ops::Sub<Scalar> for Scalar {
            type Output = Scalar;

            fn sub(self, other: Scalar) -> Scalar {
                &self - &other
            }
        }

        // **********************
        // Scalar Multiplication
        // **********************

        impl<'a, 'b> std::ops::Mul<&'b Scalar> for &'a Scalar {
            type Output = Scalar;

            fn mul(self, other: &'b Scalar) -> Scalar {
                Scalar((&self.0 * &other.0) % $p)
            }
        }

        impl<'b> std::ops::Mul<&'b Scalar> for Scalar {
            type Output = Scalar;

            fn mul(self, other: &'b Scalar) -> Scalar {
                Scalar((&self.0 * &other.0) % $p)
            }
        }

        impl<'a, 'b> std::ops::Mul<Scalar> for &'a Scalar {
            type Output = Scalar;

            fn mul(self, other: Scalar) -> Scalar {
                Scalar((&self.0 * &other.0) % $p)
            }
        }

        impl std::ops::Mul<Scalar> for Scalar {
            type Output = Scalar;

            fn mul(self, other: Scalar) -> Scalar {
                Scalar((&self.0 * &other.0) % $p)
            }
        }
    };
}
