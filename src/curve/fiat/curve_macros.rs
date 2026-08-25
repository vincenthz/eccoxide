#[doc(hidden)]
#[macro_export]
macro_rules! fiat_define_weierstrass_curve {
    ($FE:ident) => {
        const A: $FE = $FE::from_bytes_unchecked(&A_BYTES);
        const B: $FE = $FE::from_bytes_unchecked(&B_BYTES);
        const B3: $FE = $FE::from_bytes_unchecked(&B3_BYTES);
        const GX: $FE = $FE::from_bytes_unchecked(&GX_BYTES);
        const GY: $FE = $FE::from_bytes_unchecked(&GY_BYTES);

        /// The Weierstrass elliptic curve object itself
        #[derive(Debug, Clone, Copy)]
        pub struct Curve;

        impl Curve {
            /// Get the group order as an array of bytes in big endian representation
            pub fn group_order(self) -> &'static [u8] {
                &ORDER_BYTES
            }

            /// Return the generator field element in affine coordinate (X,Y)
            pub fn generator() -> (&'static $FE, &'static $FE) {
                (&GX, &GY)
            }
        }

        impl WeierstrassCurve for Curve {
            type FieldElement = $FE;

            const A: Self::FieldElement = A;
            const B: Self::FieldElement = B;
            const B3: Self::FieldElement = B3;
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! fiat_define_weierstrass_curve_am3 {
    ($FE:ident) => {
        impl WeierstrassCurveAM3 for Curve {}

        impl Point {
            fn add_or_double<'b>(&self, other: &'b Point) -> Point {
                Point(self.0.add_or_double_am3::<Curve>(&other.0))
            }
            fn scale<'b>(&self, other: &'b Scalar) -> Self {
                Point(self.0.scale_am3_ct::<Curve>(&other.to_bytes()))
            }

            /// Constant-time multiplication of the curve generator: `scalar * G`.
            ///
            /// Uses the precomputed fixed-base comb table, which is substantially
            /// faster than the general `&Point * &Scalar` path.
            pub fn mul_base(scalar: &Scalar) -> Self {
                #[cfg(feature = "table")]
                return Point(projective::Point::<$FE>::mul_base_table_am3::<_, Curve>(
                    generator_comb(),
                    &scalar.to_bytes(),
                ));
                #[cfg(not(feature = "table"))]
                return &Point::GENERATOR * scalar;
            }

            /// Variable-time scalar multiplication.
            ///
            /// Faster than the constant-time `*` operator, but its running time depends
            /// on the scalar; only use it when the scalar is public.
            pub fn mul_vartime(&self, other: &Scalar) -> Self {
                Point(self.0.scale_am3::<Curve>(&other.to_bytes()))
            }
        }

        impl $crate::curve::group::CurveGroup for Point {
            type Scalar = Scalar;

            const IDENTITY: Self = Point::INFINITY;
            const GENERATOR: Self = Point::GENERATOR;

            fn double(&self) -> Self {
                Point(self.0.double_am3::<Curve>())
            }
            fn mul_base(scalar: &Scalar) -> Self {
                Point::mul_base(scalar)
            }
            fn mul_vartime(&self, scalar: &Scalar) -> Self {
                Point::mul_vartime(self, scalar)
            }
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! fiat_define_weierstrass_curve_a0 {
    ($FE:ident) => {
        impl WeierstrassCurveA0 for Curve {}

        impl Point {
            fn add_or_double<'b>(&self, other: &'b Point) -> Point {
                Point(self.0.add_or_double_a0::<Curve>(&other.0))
            }
            fn scale<'b>(&self, other: &'b Scalar) -> Self {
                Point(self.0.scale_a0_ct::<Curve>(&other.to_bytes()))
            }

            /// Constant-time multiplication of the curve generator: `scalar * G`.
            ///
            /// Uses the precomputed fixed-base comb table, which is substantially
            /// faster than the general `&Point * &Scalar` path.
            pub fn mul_base(scalar: &Scalar) -> Self {
                #[cfg(feature = "table")]
                return Point(projective::Point::<$FE>::mul_base_table_a0::<_, Curve>(
                    generator_comb(),
                    &scalar.to_bytes(),
                ));
                #[cfg(not(feature = "table"))]
                return &Point::GENERATOR * scalar;
            }

            /// Variable-time scalar multiplication.
            ///
            /// Faster than the constant-time `*` operator, but its running time depends
            /// on the scalar; only use it when the scalar is public.
            pub fn mul_vartime(&self, other: &Scalar) -> Self {
                Point(self.0.scale_a0::<Curve>(&other.to_bytes()))
            }
        }

        impl $crate::curve::group::CurveGroup for Point {
            type Scalar = Scalar;

            const IDENTITY: Self = Point::INFINITY;
            const GENERATOR: Self = Point::GENERATOR;

            fn double(&self) -> Self {
                Point(self.0.double_a0::<Curve>())
            }
            fn mul_base(scalar: &Scalar) -> Self {
                Point::mul_base(scalar)
            }
            fn mul_vartime(&self, scalar: &Scalar) -> Self {
                Point::mul_vartime(self, scalar)
            }
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! fiat_define_weierstrass_points {
    ($FE:ident) => {
        /// Affine Point on the curve of type (X,Y)
        ///
        /// Note that this representation cannot handle the point at infinity
        #[derive(Clone, Debug, PartialEq, Eq)]
        pub struct PointAffine(affine::Point<$FE>);

        /// Point on the curve using a more optimised representation
        ///
        /// This implementation used projective coordinate (X:Y:Z)
        #[derive(Clone, Debug, PartialEq, Eq)]
        pub struct Point(projective::Point<$FE>);

        /// Fixed-base comb table for the generator, built once on first use
        /// from the statically embedded `COMB_TABLE` constant.
        #[cfg(feature = "table")]
        fn generator_comb() -> &'static [[projective::Point<$FE>; 16]; COMB_WINDOWS] {
            // The comb table is held boxed (see `build_comb_table`) so it is
            // never materialized on the stack; deref to a plain array reference.
            static V: std::sync::OnceLock<Box<[[projective::Point<$FE>; 16]; COMB_WINDOWS]>> =
                std::sync::OnceLock::new();
            &**V.get_or_init(|| {
                projective::Point::<$FE>::build_comb_table(&COMB_TABLE, $FE::from_bytes_unchecked)
            })
        }

        impl PointAffine {
            /// Curve generator point in affine coordinate
            pub const GENERATOR: Self = PointAffine(affine::Point { x: GX, y: GY });

            /// Try to create an affine point with X, Y coordinates.
            ///
            /// check if the equation y^2 = x^3 + a*x + b (mod p) holds for this curve, if it doesn't
            /// None is returned
            pub fn from_coordinate(x: &FieldElement, y: &FieldElement) -> Option<Self> {
                affine::Point::from_coordinate::<Curve>(x, y).map(PointAffine)
            }

            /// Return the tuple of coordinate (x, y) associated with this
            /// affine point
            pub const fn to_coordinate(&self) -> (&FieldElement, &FieldElement) {
                (&self.0.x, &self.0.y)
            }

            /// Double the affine point Self
            ///
            /// This is equivalent to Self + Self at the mathematic level,
            /// but is implemented more quickly than the normal addition
            /// of double possibly arbitrary point
            pub fn double(&self) -> PointAffine {
                PointAffine(affine::Point::double::<Curve>(&self.0))
            }

            /// Turn an affine point into the X component and the sign of the Y component
            ///
            /// This is often refered as point compression, and related to the fact there
            /// two point on the curve for a valid x component as (x,y) and (x,-y), unless
            /// y is 0. So it is sufficient to know just the sign of y to know which point
            /// is in use for a given x component
            pub fn compress(&self) -> (&FieldElement, Sign) {
                self.0.compress()
            }

            /// Try to create an affine point given a X component and the sign
            /// of the Y component.
            ///
            /// This is often refered as point decompression.
            ///
            /// It is done in constant time and does not branch on the recovered point value
            pub fn decompress(x: &FieldElement, sign: Sign) -> $crate::mp::ct::CtOption<Self> {
                affine::Point::decompress::<Curve>(x, sign).map(PointAffine)
            }
        }

        impl<'a, 'b> core::ops::Add<&'b PointAffine> for &'a PointAffine {
            type Output = PointAffine;
            fn add(self, other: &'b PointAffine) -> PointAffine {
                PointAffine(&self.0 + &other.0)
            }
        }

        impl Point {
            /// Curve generator point
            pub const GENERATOR: Self = Point(projective::Point {
                x: GX,
                y: GY,
                z: FieldElement::ONE,
            });

            /// Point at infinity, used as additive zero
            pub const INFINITY: Self = Point(projective::Point::<$FE>::INFINITY);

            /// Convert an affine point to optimised point representation
            ///
            /// In projective coordinate it means, (X,Y) => (X:Y:1)
            pub fn from_affine(p: &PointAffine) -> Self {
                Point(projective::Point::from_affine(&p.0))
            }

            /// Convert a point to the affine point
            ///
            /// In projective coordinate it means, (X:Y:Z) => (X/Z, Y/Z)
            pub fn to_affine(&self) -> Option<PointAffine> {
                self.0.to_affine().map(PointAffine)
            }

            /// Constant-time variant of [`Self::to_affine`]: returns a
            /// `CtOption` marked not-present for the point at infinity,
            /// without branching on the point value
            pub fn to_affine_ct(&self) -> $crate::mp::ct::CtOption<PointAffine> {
                let (present, p) = self.0.to_affine_ct().into_parts();
                $crate::mp::ct::CtOption::from((present, PointAffine(p)))
            }

            /// Constant-time affine x-coordinate: like [`Self::to_affine_ct`]
            /// but computes only `x/z`, skipping the y-coordinate. Marked
            /// not-present for the point at infinity.
            pub fn to_affine_x_ct(&self) -> $crate::mp::ct::CtOption<FieldElement> {
                self.0.to_affine_x_ct()
            }

            /// Normalize the point, keeping the same representation
            ///
            /// In projective coordinate it means, (X:Y:Z) => (X/Z:Y/Z:1)
            pub fn normalize(&mut self) {
                self.0.normalize()
            }
        }

        impl From<PointAffine> for Point {
            fn from(p: PointAffine) -> Self {
                Point(projective::Point::from_affine(&p.0))
            }
        }

        impl From<&PointAffine> for Point {
            fn from(p: &PointAffine) -> Self {
                Point(projective::Point::from_affine(&p.0))
            }
        }

        // *************
        // Point Negation
        // *************

        impl core::ops::Neg for Point {
            type Output = Point;

            fn neg(self) -> Self::Output {
                Point(self.0.neg())
            }
        }

        impl<'a> core::ops::Neg for &'a Point {
            type Output = Point;

            fn neg(self) -> Self::Output {
                Point(self.0.clone().neg())
            }
        }

        // *************
        // Point Scaling
        // *************

        // note that scalar multiplication is really defined for arbitrary scalar
        // (of any size), not just the *field element* scalar defined in F(p).
        // this semantic abuse makes it easier to use.

        impl<'a, 'b> core::ops::Mul<&'b Scalar> for &'a Point {
            type Output = Point;

            fn mul(self, other: &'b Scalar) -> Point {
                self.scale(other)
            }
        }

        impl<'b> core::ops::Mul<&'b Scalar> for Point {
            type Output = Point;

            fn mul(self, other: &'b Scalar) -> Point {
                self.scale(other)
            }
        }

        impl<'a, 'b> core::ops::Mul<&'b Point> for &'a Scalar {
            type Output = Point;

            fn mul(self, other: &'b Point) -> Point {
                other * self
            }
        }

        // **************
        // Point Addition
        // **************

        impl<'a, 'b> core::ops::Add<&'b Point> for &'a Point {
            type Output = Point;

            fn add(self, other: &'b Point) -> Point {
                self.add_or_double(other)
            }
        }

        impl<'b> core::ops::Add<&'b Point> for Point {
            type Output = Point;

            fn add(self, other: &'b Point) -> Point {
                &self + other
            }
        }

        impl<'a> core::ops::Add<Point> for &'a Point {
            type Output = Point;

            fn add(self, other: Point) -> Point {
                self + &other
            }
        }

        impl core::ops::Add<Point> for Point {
            type Output = Point;

            fn add(self, other: Point) -> Point {
                &self + &other
            }
        }

        impl<'a, 'b> core::ops::Sub<&'b Point> for &'a Point {
            type Output = Point;

            fn sub(self, other: &'b Point) -> Point {
                self + (-other)
            }
        }

        impl<'b> core::ops::Sub<&'b Point> for Point {
            type Output = Point;

            fn sub(self, other: &'b Point) -> Point {
                &self - other
            }
        }

        impl<'a> core::ops::Sub<Point> for &'a Point {
            type Output = Point;

            fn sub(self, other: Point) -> Point {
                self - &other
            }
        }

        impl core::ops::Sub<Point> for Point {
            type Output = Point;

            fn sub(self, other: Point) -> Point {
                &self - &other
            }
        }
    };
}

/// Emit the point compression / decompression unit tests for a fiat curve.
///
/// Just like [`fiat_field_unittest`](crate::fiat_field_unittest) factors the
/// field-element tests, this factors the affine-point tests so every curve
/// gets the same coverage. Invoke it from inside a `mod point { ... }` of the
/// curve's `tests` module; it pulls `PointAffine`/`FieldElement` from the
/// enclosing curve module and `Sign` from the crate.
#[doc(hidden)]
#[macro_export]
macro_rules! fiat_curve_point_unittest {
    () => {
        use super::super::{FieldElement, PointAffine};
        use $crate::curve::field::Sign;

        /// `compress` followed by `decompress` recovers the original point for
        /// both possible y-signs; requesting the opposite sign yields the
        /// negated point.
        #[test]
        fn compress_decompress_roundtrip() {
            for p in [PointAffine::GENERATOR, PointAffine::GENERATOR.double()] {
                let (x, sign) = p.compress();
                let x = x.clone();

                let recovered = PointAffine::decompress(&x, sign)
                    .into_option()
                    .expect("valid x must decompress");
                assert_eq!(recovered, p);

                let other = if sign == Sign::Positive {
                    Sign::Negative
                } else {
                    Sign::Positive
                };
                let neg = PointAffine::decompress(&x, other)
                    .into_option()
                    .expect("valid x must decompress");
                assert_eq!(neg.to_coordinate().1, &(-recovered.to_coordinate().1));
            }
        }

        /// An `x` that is not a valid compressed x-coordinate must decompress
        /// to a not-present `CtOption`.
        #[test]
        fn decompress_rejects_invalid_x() {
            // find the first small x that is not on the curve
            let mut i = 0u64;
            let x = loop {
                let x = FieldElement::from_u64(i);
                if PointAffine::decompress(&x, Sign::Positive)
                    .into_option()
                    .is_none()
                {
                    break x;
                }
                i += 1;
            };
            assert!(PointAffine::decompress(&x, Sign::Positive)
                .into_option()
                .is_none());
        }
    };
}
/*

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bytes() {
        let b: [u8; 32] = [
            1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 7, 8, 9, 0, 1,
            2, 3, 4,
        ];
        let s = Scalar::from_bytes(&b).unwrap();
        assert_eq!(b, s.to_bytes())
    }
    #[test]
    fn bytes_u64() {
        let b: [u8; 32] = [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0xce,
        ];

        let s1 = Scalar::from_bytes(&b).unwrap();
        let s2 = Scalar::from_u64(0xce);
        assert_eq!(s1, s2)
    }
}

*/
