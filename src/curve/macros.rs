//!
//! macros to generate various types (Scalar, PointAffine, Point) given different curve properties.
//!
//! Reference reading:
//!
//! * [Complete addition formulas for prime order elliptic curves](https://eprint.iacr.org/2015/1060.pdf) (1)
//! * Handbook of Elliptic and Hyperelliptic Curve Cryptography - Chapter 13
//! * [NIST.SP.800-186](https://csrc.nist.gov/publications/detail/sp/800-186/draft) : Appendix D & E

#[doc(hidden)]
#[macro_export]
macro_rules! scalar_impl {
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

#[doc(hidden)]
#[macro_export]
macro_rules! point_impl {
    ($gx: expr, $gy: expr) => {
        /// Affine Point on the curve
        #[derive(Clone, Debug, PartialEq, Eq)]
        pub struct PointAffine {
            x: FieldElement,
            y: FieldElement,
        }

        /// Point on the curve
        #[derive(Clone, Debug)]
        pub struct Point {
            x: FieldElement,
            y: FieldElement,
            z: FieldElement,
        }

        impl PartialEq for Point {
            fn eq(&self, other: &Point) -> bool {
                self.to_affine() == other.to_affine()
            }
        }

        impl Eq for Point {}

        lazy_static! {
            static ref A: FieldElement = FieldElement(BigUint::from_bytes_be(&A_BYTES));
            static ref B: FieldElement = FieldElement(BigUint::from_bytes_be(&B_BYTES));
            static ref B3: FieldElement =
                FieldElement(BigUint::from_bytes_be(&B_BYTES) * BigUint::from_bytes_be(&[3]));
            static ref GX: FieldElement = FieldElement(BigUint::from_bytes_be(&GX_BYTES));
            static ref GY: FieldElement = FieldElement(BigUint::from_bytes_be(&GY_BYTES));
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
            pub fn from_coordinate(x: &FieldElement, y: &FieldElement) -> Option<Self> {
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

            pub fn to_coordinate(&self) -> (&FieldElement, &FieldElement) {
                (&self.x, &self.y)
            }

            pub fn double(&self) -> PointAffine {
                let PointAffine {
                    x: ref x1,
                    y: ref y1,
                } = self;
                let l = (FieldElement::from_u64(3) * (x1 * x1) + &*A)
                    * (FieldElement::from_u64(2) * y1).inverse().unwrap();
                let l2 = &l * &l;
                let x3 = l2 - FieldElement::from_u64(2) * x1;
                let y3 = l * (x1 - &x3) - y1;
                PointAffine { x: x3, y: y3 }
            }

            pub fn compress(&self) -> (&FieldElement, bool) {
                (&self.x, self.y.high_bit_set())
            }

            pub fn decompress(x: &FieldElement, bit: bool) -> Option<Self> {
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
                    z: FieldElement::one(),
                }
            }

            /// Point at infinity
            pub fn infinity() -> Self {
                Point {
                    x: FieldElement::zero(),
                    y: FieldElement::one(),
                    z: FieldElement::zero(),
                }
            }

            pub fn from_affine(p: &PointAffine) -> Self {
                Point {
                    x: p.x.clone(),
                    y: p.y.clone(),
                    z: FieldElement::one(),
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
                self.z = FieldElement::one()
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
                    z: FieldElement::one(),
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

        impl<'a, 'b> std::ops::Mul<&'b Scalar> for &'a Point {
            type Output = Point;

            fn mul(self, other: &'b Scalar) -> Point {
                self.scalar_mul_daa_limbs8(&other.to_bytes())
            }
        }

        impl<'a, 'b> std::ops::Mul<&'b Point> for &'a Scalar {
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
    () => {
        #[test]
        fn point_add_infinity() {
            let p = &Point::generator() * &Scalar::from_u64(1245);
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
            let mut scalar = [0u8; Scalar::SIZE_BYTES];
            scalar[Scalar::SIZE_BYTES - 1] = 0x2;
            scalar[Scalar::SIZE_BYTES - 10] = 0xd6;
            let s = Scalar::from_bytes(&scalar).unwrap();
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
            let p2got = &p1 * &Scalar::from_u64(2);
            let p4got = &p1 * &Scalar::from_u64(4);

            let p6got = &p1 * &Scalar::from_u64(6);
            let p8got = &p1 * &Scalar::from_u64(8);

            assert_eq!(p2, p2got);
            assert_eq!(p4, p4got);
            assert_eq!(p6, p6got);
            assert_eq!(p8, p8got);
            assert_eq!(&p2got * &Scalar::from_u64(4), p8got);

            for b in &[4u64, 8, 10, 11, 39] {
                let g = &Point::generator() * &Scalar::from_u64(*b);
                for p in &[34u64, 56, 791, 12492124] {
                    let p1: u64 = p / 2;
                    let p2: u64 = p - p1;

                    let r = &g * &Scalar::from_u64(*p);
                    let r1 = &g * &Scalar::from_u64(p1);
                    let r2 = &g * &Scalar::from_u64(p2);
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
