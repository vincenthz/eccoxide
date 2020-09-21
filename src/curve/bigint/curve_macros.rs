//!
//! macros to generate various types ($SCALAR, PointAffine, Point) given different curve properties.
//!
//! Reference reading:
//!
//! * [Complete addition formulas for prime order elliptic curves](https://eprint.iacr.org/2015/1060.pdf) (1)
//! * Handbook of Elliptic and Hyperelliptic Curve Cryptography - Chapter 13
//! * [NIST.SP.800-186](https://csrc.nist.gov/publications/detail/sp/800-186/draft) : Appendix D & E

#[doc(hidden)]
#[macro_export]
macro_rules! point_impl {
    ($FE:ident, $SCALAR:ident, $gx:expr, $gy:expr) => {
        /// Affine Point on the curve
        #[derive(Clone, Debug, PartialEq, Eq)]
        pub struct PointAffine {
            x: $FE,
            y: $FE,
        }

        /// Point on the curve
        #[derive(Clone, Debug)]
        pub struct Point {
            x: $FE,
            y: $FE,
            z: $FE,
        }

        impl PartialEq for Point {
            fn eq(&self, other: &Point) -> bool {
                self.to_affine() == other.to_affine()
            }
        }

        impl Eq for Point {}

        lazy_static! {
            static ref A: $FE = $FE(BigUint::from_bytes_be(&A_BYTES));
            static ref B: $FE = $FE(BigUint::from_bytes_be(&B_BYTES));
            static ref B3: $FE =
                $FE(BigUint::from_bytes_be(&B_BYTES) * BigUint::from_bytes_be(&[3]));
            static ref GX: $FE = $FE(BigUint::from_bytes_be(&GX_BYTES));
            static ref GY: $FE = $FE(BigUint::from_bytes_be(&GY_BYTES));
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
            pub fn from_coordinate(x: &$FE, y: &$FE) -> Option<Self> {
                let y2 = y * y;
                let x3 = x.power(3);
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

            pub fn to_coordinate(&self) -> (&$FE, &$FE) {
                (&self.x, &self.y)
            }

            pub fn double(&self) -> PointAffine {
                let PointAffine {
                    x: ref x1,
                    y: ref y1,
                } = self;
                let l = ($FE::from_u64(3) * (x1 * x1) + &*A)
                    * ($FE::from_u64(2) * y1).inverse().unwrap();
                let l2 = &l * &l;
                let x3 = l2 - $FE::from_u64(2) * x1;
                let y3 = l * (x1 - &x3) - y1;
                PointAffine { x: x3, y: y3 }
            }

            pub fn compress(&self) -> (&$FE, bool) {
                (&self.x, self.y.high_bit_set())
            }

            pub fn decompress(x: &$FE, bit: bool) -> Option<Self> {
                // Y^2 = X^3 - A*X + b
                let yy = x.power(3) + (&*A * x) + &*B;
                let y = yy.sqrt()?;
                let x = x.clone();
                if bit == y.high_bit_set() {
                    Some(PointAffine { x, y })
                } else {
                    Some(PointAffine { x, y: -y })
                }
            }
        }

        impl<'a, 'b> std::ops::Add<&'b PointAffine> for &'a PointAffine {
            type Output = PointAffine;
            fn add(self, other: &'b PointAffine) -> PointAffine {
                let PointAffine {
                    x: ref x1,
                    y: ref y1,
                } = &self;
                let PointAffine {
                    x: ref x2,
                    y: ref y2,
                } = &other;
                let l = (y1 - y2) * (x1 - x2).inverse().expect("inverse exist");
                let l2 = &l * &l;
                let x3 = l2 - x1 - x2;
                let y3 = &l * (x1 - &x3) - y1;
                PointAffine { x: x3, y: y3 }
            }
        }

        impl Point {
            /// Curve generator point
            pub fn generator() -> Self {
                Point {
                    x: GX.clone(),
                    y: GY.clone(),
                    z: $FE::one(),
                }
            }

            /// Point at infinity
            pub fn infinity() -> Self {
                Point {
                    x: $FE::zero(),
                    y: $FE::one(),
                    z: $FE::zero(),
                }
            }

            pub fn from_affine(p: &PointAffine) -> Self {
                Point {
                    x: p.x.clone(),
                    y: p.y.clone(),
                    z: $FE::one(),
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

            pub fn normalize(&mut self) {
                let zinv = self.z.inverse().unwrap();

                self.x = &self.x * &zinv;
                self.y = &self.y * &zinv;
                self.z = $FE::one()
            }

            fn add_different<'b>(&self, other: &'b Point) -> Point {
                assert!(self != other);

                let Point {
                    x: ref x1,
                    y: ref y1,
                    z: ref z1,
                } = &self;
                let Point {
                    x: ref x2,
                    y: ref y2,
                    z: ref z2,
                } = &other;

                // Algorithm 1 from (1) - addition formula for arbitrary a
                //
                // ```magma
                // ADD := function ( X1 , Y1 , Z1 , X2 , Y2 , Z2 ,a , b3 )
                //     t0 := X1 * X2 ; t1 := Y1 * Y2 ; t2 := Z1 * Z2 ;
                //     t3 := X1 + Y1 ; t4 := X2 + Y2 ; t3 := t3 * t4 ;
                //     t4 := t0 + t1 ; t3 := t3 - t4 ; t4 := X1 + Z1 ;
                //     t5 := X2 + Z2 ; t4 := t4 * t5 ; t5 := t0 + t2 ;
                //     t4 := t4 - t5 ; t5 := Y1 + Z1 ; X3 := Y2 + Z2 ;
                //     t5 := t5 * X3 ; X3 := t1 + t2 ; t5 := t5 - X3 ;
                //     Z3 := a * t4 ; X3 := b3 * t2 ; Z3 := X3 + Z3 ;
                //     X3 := t1 - Z3 ; Z3 := t1 + Z3 ; Y3 := X3 * Z3 ;
                //     t1 := t0 + t0 ; t1 := t1 + t0 ; t2 := a * t2 ;
                //     t4 := b3 * t4 ; t1 := t1 + t2 ; t2 := t0 - t2 ;
                //     t2 := a * t2 ; t4 := t4 + t2 ; t0 := t1 * t4 ;
                //     Y3 := Y3 + t0 ; t0 := t5 * t4 ; X3 := t3 * X3 ;
                //     X3 := X3 - t0 ; t0 := t3 * t1 ; Z3 := t5 * Z3 ;
                //     Z3 := Z3 + t0 ;
                //     return X3 , Y3 , Z3 ;
                // end function ;
                // ```

                let t0 = x1 * x2;
                let t1 = y1 * y2;
                let t2 = z1 * z2;
                let t3 = x1 + y1;
                let t4 = x2 + y2;
                let t3 = t3 * t4;
                let t4 = &t0 + &t1;
                let t3 = t3 - &t4;
                let t4 = x1 + z1;
                let t5 = x2 + z2;
                let t4 = t4 * &t5;
                let t5 = &t0 + &t2;
                let t4 = t4 - &t5;
                let t5 = y1 + z1;
                let x3 = y2 + z2;
                let t5 = t5 * &x3;
                let x3 = &t1 + &t2;
                let t5 = t5 - &x3;
                let z3 = &*A * &t4;
                let x3 = &*B3 * &t2;
                let z3 = &x3 + z3;
                let x3 = &t1 - &z3;
                let z3 = &t1 + &z3;
                let y3 = &x3 * &z3;
                let t1 = t0.double();
                let t1 = t1 + &t0;
                let t2 = &*A * t2;
                let t4 = &*B3 * &t4;
                let t1 = t1 + &t2;
                let t2 = &t0 - t2;
                let t2 = &*A * t2;
                let t4 = t4 + &t2;
                let t0 = &t1 * &t4;
                let y3 = y3 + t0;
                let t0 = &t5 * &t4;
                let x3 = &t3 * x3;
                let x3 = x3 - &t0;
                let t0 = &t3 * &t1;
                let z3 = &t5 * &z3;
                let z3 = z3 + t0;

                Point {
                    x: x3,
                    y: y3,
                    z: z3,
                }
            }

            pub fn double(&self) -> Self {
                let Point {
                    x: ref x,
                    y: ref y,
                    z: ref z,
                } = &self;

                // Algorithm 3 from (1) - doubling formula for arbitrary a
                // ```magma
                // DBL := function (X ,Y ,Z ,a , b3 )
                //    t0 := X ^2; t1 := Y ^2; t2 := Z ^2;
                //    t3 := X * Y ; t3 := t3 + t3 ; Z3 := X * Z ;
                //    Z3 := Z3 + Z3 ; X3 := a * Z3 ; Y3 := b3 * t2 ;
                //    Y3 := X3 + Y3 ; X3 := t1 - Y3 ; Y3 := t1 + Y3 ;
                //    Y3 := X3 * Y3 ; X3 := t3 * X3 ; Z3 := b3 * Z3 ;
                //    t2 := a * t2 ; t3 := t0 - t2 ; t3 := a * t3 ;
                //    t3 := t3 + Z3 ; Z3 := t0 + t0 ; t0 := Z3 + t0 ;
                //    t0 := t0 + t2 ; t0 := t0 * t3 ; Y3 := Y3 + t0 ;
                //    t2 := Y * Z ; t2 := t2 + t2 ; t0 := t2 * t3 ;
                //    X3 := X3 - t0 ; Z3 := t2 * t1 ; Z3 := Z3 + Z3 ;
                //    Z3 := Z3 + Z3 ;
                //    return X3 , Y3 , Z3 ;
                // end function ;
                // ```
                let t0 = x * x;
                let t1 = y * y;
                let t2 = z * z;
                let t3 = x * y;
                let t3 = t3.double();
                let z3 = x * z;
                let z3 = &z3 + &z3;
                let x3 = &*A * &z3;
                let y3 = &*B3 * &t2;
                let y3 = &x3 + &y3;
                let x3 = &t1 - &y3;
                let y3 = &t1 + &y3;
                let y3 = &x3 * &y3;
                let x3 = &t3 * &x3;
                let z3 = &*B3 * &z3;
                let t2 = &*A * &t2;
                let t3 = &t0 - &t2;
                let t3 = &*A * &t3;
                let t3 = &t3 + &z3;
                let z3 = &t0 + &t0;
                let t0 = &z3 + &t0;
                let t0 = &t0 + &t2;
                let t0 = &t0 * &t3;
                let y3 = &y3 + &t0;
                let t2 = y * z;
                let t2 = &t2 + &t2;
                let t0 = &t2 * &t3;
                let x3 = &x3 - &t0;
                let z3 = &t2 * &t1;
                let z3 = &z3 + &z3;
                let z3 = &z3 + &z3;

                Point {
                    x: x3,
                    y: y3,
                    z: z3,
                }
            }

            fn add_or_double<'b>(&self, other: &'b Point) -> Point {
                if self == other {
                    self.double()
                } else {
                    self.add_different(other)
                }
            }

            /// scalar multiplication : `n * self` with double-and-add algorithm with increasing index
            fn scalar_mul_daa_limbs8(&self, n: &[u8]) -> Self {
                let mut a = self.clone();
                let mut q = Point::infinity();

                for digit in n.iter().rev() {
                    for i in 0..8 {
                        if digit & (1 << i) != 0 {
                            q = &q + &a;
                        }
                        a = a.double()
                    }
                }
                q
            }
        }

        impl From<PointAffine> for Point {
            fn from(p: PointAffine) -> Self {
                Point {
                    x: p.x,
                    y: p.y,
                    z: $FE::one(),
                }
            }
        }

        impl From<&PointAffine> for Point {
            fn from(p: &PointAffine) -> Self {
                Point::from_affine(p)
            }
        }

        // *************
        // Point Negation
        // *************

        impl std::ops::Neg for Point {
            type Output = Point;

            fn neg(self) -> Self::Output {
                Point {
                    x: self.x,
                    y: -self.y,
                    z: self.z,
                }
            }
        }

        impl<'a> std::ops::Neg for &'a Point {
            type Output = Point;

            fn neg(self) -> Self::Output {
                Point {
                    x: self.x.clone(),
                    y: -&self.y,
                    z: self.z.clone(),
                }
            }
        }

        // *************
        // Point Scaling
        // *************

        // note that scalar multiplication is really defined for arbitrary scalar
        // (of any size), not just the *field element* scalar defined in F(p).
        // this semantic abuse makes it easier to use.

        impl<'a, 'b> std::ops::Mul<&'b $SCALAR> for &'a Point {
            type Output = Point;

            fn mul(self, other: &'b $SCALAR) -> Point {
                self.scalar_mul_daa_limbs8(&other.to_bytes())
            }
        }

        impl<'a, 'b> std::ops::Mul<&'b Point> for &'a $SCALAR {
            type Output = Point;

            fn mul(self, other: &'b Point) -> Point {
                other * self
            }
        }

        // **************
        // Point Addition
        // **************

        impl<'a, 'b> std::ops::Add<&'b Point> for &'a Point {
            type Output = Point;

            fn add(self, other: &'b Point) -> Point {
                self.add_or_double(other)
            }
        }

        impl<'b> std::ops::Add<&'b Point> for Point {
            type Output = Point;

            fn add(self, other: &'b Point) -> Point {
                &self + other
            }
        }

        impl<'a> std::ops::Add<Point> for &'a Point {
            type Output = Point;

            fn add(self, other: Point) -> Point {
                self + &other
            }
        }

        impl std::ops::Add<Point> for Point {
            type Output = Point;

            fn add(self, other: Point) -> Point {
                &self + &other
            }
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! test_scalar_arithmetic {
    ($scalar: ident) => {
        use super::$scalar;

        #[test]
        fn negate() {
            assert_eq!($scalar::one() + (-$scalar::one()), $scalar::zero())
        }

        #[test]
        fn high() {
            assert!(!($scalar::one().high_bit_set()), "1");
            assert!((-$scalar::one()).high_bit_set(), "-1");
        }

        #[test]
        fn inverse() {
            assert_eq!(
                $scalar::one() * $scalar::one().inverse().unwrap(),
                $scalar::one()
            );

            let mut v = $scalar::one() + $scalar::one();
            for _ in 0..100 {
                assert_eq!(&v * v.inverse().unwrap(), $scalar::one());
                v = v + $scalar::one();
            }

            for i in 1..16 {
                let s = $scalar::from_u64(i * 1048);
                let sinv = s.inverse().unwrap();
                let r = &s * &sinv;
                assert_eq!(r, $scalar::one());
            }
        }

        #[test]
        fn sqrt() {
            let y = $scalar::one().sqrt().unwrap();
            assert_eq!(&y * &y, $scalar::one());
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! test_point_arithmetic {
    ($SCALAR:ident) => {
        #[test]
        fn point_add_infinity() {
            let p = &Point::generator() * &$SCALAR::from_u64(1245);
            assert_eq!(Point::generator() + Point::infinity(), Point::generator());
            assert_eq!(&p + Point::infinity(), p);
            assert_eq!(Point::infinity() + &p, p);
        }

        #[test]
        fn point_affine_projective() {
            assert_eq!(
                Point::from(Point::generator().to_affine().unwrap()),
                Point::generator()
            )
        }

        #[test]
        fn point_mul_and_inv() {
            let mut scalar = [0u8; $SCALAR::SIZE_BYTES];
            scalar[$SCALAR::SIZE_BYTES - 1] = 0x2;
            scalar[$SCALAR::SIZE_BYTES - 10] = 0xd6;
            let s = $SCALAR::from_bytes(&scalar).unwrap();
            let sinv = s.inverse().unwrap();
            let p = &Point::generator() * &s;
            let p2 = &p * &sinv;
            assert_eq!(&p2, &Point::generator())
        }

        #[test]
        fn point_mul() {
            let p1 = Point::generator();

            let p2 = p1.double();
            let p4 = p2.double();
            let p6 = &p2 + &p4;
            let p8 = p4.double();
            let p2got = &p1 * &$SCALAR::from_u64(2);
            let p4got = &p1 * &$SCALAR::from_u64(4);

            let p6got = &p1 * &$SCALAR::from_u64(6);
            let p8got = &p1 * &$SCALAR::from_u64(8);

            assert_eq!(p2, p2got);
            assert_eq!(p4, p4got);
            assert_eq!(p6, p6got);
            assert_eq!(p8, p8got);
            assert_eq!(&p2got * &$SCALAR::from_u64(4), p8got);

            for b in &[4u64, 8, 10, 11, 39] {
                let g = &Point::generator() * &$SCALAR::from_u64(*b);
                for p in &[34u64, 56, 791, 12492124] {
                    let p1: u64 = p / 2;
                    let p2: u64 = p - p1;

                    let r = &g * &$SCALAR::from_u64(*p);
                    let r1 = &g * &$SCALAR::from_u64(p1);
                    let r2 = &g * &$SCALAR::from_u64(p2);
                    let rprim = r1 + r2;

                    assert_eq!(r, rprim, "(p1 + p2) X == p1 X + p2 X");
                }
            }
        }

        #[test]
        fn point_serialization() {
            let p = PointAffine::generator();
            let (x, ysign) = p.compress();
            assert_eq!(p, PointAffine::decompress(x, ysign).unwrap());
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! bigint_prime_curve {
    ($m: ident, $szfe: expr) => {
        pub mod $m {
            use crate::curve::bigint::maths::{mod_inverse, tonelli_shanks};
            use crate::params::sec2::$m::*;
            use crate::{bigint_scalar_impl, point_impl};
            use lazy_static;
            use num_bigint::BigUint;
            use num_traits::{cast::ToPrimitive, identities::One};

            lazy_static! {
                static ref P: BigUint = BigUint::from_bytes_be(&P_BYTES);
                static ref PMOD4: u32 = {
                    let pmodded = &*P & BigUint::from(0b11u32);
                    pmodded.to_u32().unwrap()
                };

                // "constant" (P + 1) / 4
                static ref PP1D4: BigUint = (&*P + BigUint::one()) / BigUint::from(4u32);

                static ref ORDER: BigUint = BigUint::from_bytes_be(&ORDER_BYTES);
                static ref OMOD4: u32 = {
                    let pmodded = &*ORDER & BigUint::from(0b11u32);
                    pmodded.to_u32().unwrap()
                };

                // "constant" (ORDER + 1) / 4
                static ref OP1D4: BigUint = (&*P + BigUint::one()) / BigUint::from(4u32);
            }
            bigint_scalar_impl!(FieldElement, &*P, $szfe, PMOD4, PP1D4);
            bigint_scalar_impl!(Scalar, &*ORDER, $szfe, OMOD4, OP1D4);
            point_impl!(FieldElement, Scalar, &*GX, &*GY);

            #[cfg(test)]
            mod tests {
                use super::*;
                use crate::{test_point_arithmetic, test_scalar_arithmetic};

                mod scalar {
                    use super::*;
                    test_scalar_arithmetic!(Scalar);
                }
                mod field_element {
                    use super::*;
                    test_scalar_arithmetic!(FieldElement);
                }
                test_point_arithmetic!(Scalar);
            }
            #[cfg(test)]
            mod bench {
                // placeholder
            }
        }
    };
}
