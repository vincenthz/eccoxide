#[doc(hidden)]
#[macro_export]
macro_rules! prime_weirstrass_curve {

macro_rules! scalar_impl {
    ($sz) => {
        #[derive(Debug, Clone)]
        pub struct Scalar([u64; $sz]);

        impl PartialEq for Scalar {
            fn eq(&self, other: &Self) -> bool {
                todo!()
            }
        }

        impl Scalar {
            /// the zero constant (additive identity)
            pub const fn zero() -> Self {
                Scalar([0u64; $sz])
            }

            /// The one constant (multiplicative identity)
            pub const fn one() -> Self {
                let mut b = [0u64; $sz];
                b[$sz - 1] = 0x1;
                Scalar(b)
            }

            /// Add another scalar
            pub fn add(&self, other: &Self) -> Self {
                todo!()
            }

            /// Self add another Scalar
            pub fn add_assign(&mut self, other: &Self) {
                todo!()
            }

            /// Negate the scalar
            pub fn negate(&mut self) {
                todo!()
            }

            /// Multiplicate 2 Scalar
            pub fn mul(&self, other: &Self) -> Self {
                todo!()
            }

            /// Get the multiplicative inverse
            ///
            /// Note that 0 doesn't have a multiplicative inverse
            pub fn inverse(&self) -> Option<Self> {
                todo!()
            }

            /*
            pub fn from_bytes(&self, ) -> Option<Self> {

            }
            pub fn to_bytes(&self) -> [u8; $sz] {
            }
            */
        }
    };
}
