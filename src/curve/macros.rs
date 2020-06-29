#[doc(hidden)]
#[macro_export]
macro_rules! scalar_impl {
    ($p: expr, $sz: expr) => {
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

#[doc(hidden)]
#[macro_export]
macro_rules! point_impl {
    ($gx: expr, $gy: expr) => {
        /// Affine Point on the curve
        #[derive(Clone, Debug)]
        pub struct PointAffine {
            x: Scalar,
            y: Scalar,
        }

        /// Point on the curve
        #[derive(Clone, Debug)]
        pub struct Point {
            x: Scalar,
            y: Scalar,
            z: Scalar,
        }

        lazy_static! {
            static ref A: Scalar = Scalar(BigUint::from_bytes_be(&A_BYTES));
            static ref B: Scalar = Scalar(BigUint::from_bytes_be(&B_BYTES));
            static ref GX: Scalar = Scalar(BigUint::from_bytes_be(&GX_BYTES));
            static ref GY: Scalar = Scalar(BigUint::from_bytes_be(&GY_BYTES));
        }

        impl PointAffine {
            /// Curve generator point
            pub fn generator() -> Self {
                PointAffine {
                    x: GX.clone(),
                    y: GY.clone(),
                }
            }

            // check if y^2 = x^3 + a*x + b (mod p) holds
            pub fn from_coordinate(x: &Scalar, y: &Scalar) -> Option<Self> {
                let y2 = y * y;
                let x3 = x * x * x;
                let ax = &*A * x;

                if y2 == x3 + ax + &*B {
                    Some(PointAffine {
                        x: x.clone(),
                        y: y.clone(),
                    })
                } else {
                    None
                }
            }

            pub fn to_coordinate(&self) -> (&Scalar, &Scalar) {
                (&self.x, &self.y)
            }
        }

        impl Point {
            /// Curve generator point
            pub fn generator() -> Self {
                Point {
                    x: GX.clone(),
                    y: GY.clone(),
                    z: Scalar::one(),
                }
            }

            /// Point at infinity
            pub fn infinity() -> Self {
                Point {
                    x: Scalar::zero(),
                    y: Scalar::one(),
                    z: Scalar::zero(),
                }
            }

            pub fn from_affine(p: &PointAffine) -> Self {
                Point {
                    x: p.x.clone(),
                    y: p.y.clone(),
                    z: Scalar::one(),
                }
            }

            pub fn to_affine(&self) -> Option<PointAffine> {
                match self.z.inverse() {
                    None => None,
                    Some(inv) => Some(PointAffine {
                        x: &self.x * &inv,
                        y: &self.y * &inv,
                    }),
                }
            }
        }

        impl From<PointAffine> for Point {
            fn from(p: PointAffine) -> Self {
                Point {
                    x: p.x,
                    y: p.y,
                    z: Scalar::one(),
                }
            }
        }

        impl From<&PointAffine> for Point {
            fn from(p: &PointAffine) -> Self {
                Point::from_affine(p)
            }
        }

        // *************
        // Point Scaling
        // *************

        impl<'a, 'b> std::ops::Mul<&'b Scalar> for &'a Point {
            type Output = Point;

            fn mul(self, _other: &'b Scalar) -> Point {
                //
                todo!()
            }
        }

        // **************
        // Point Addition
        // **************

        impl<'a, 'b> std::ops::Add<&'b Point> for &'a Point {
            type Output = Point;

            fn add(self, _other: &'b Point) -> Point {
                //
                todo!()
            }
        }
    };
}
