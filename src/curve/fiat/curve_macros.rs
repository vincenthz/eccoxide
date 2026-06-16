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
            /// This is often refered as point decompression
            pub fn decompress(x: &FieldElement, sign: Sign) -> Option<Self> {
                affine::Point::decompress::<Curve>(x, sign).map(PointAffine)
            }
        }

        impl<'a, 'b> std::ops::Add<&'b PointAffine> for &'a PointAffine {
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

        impl std::ops::Neg for Point {
            type Output = Point;

            fn neg(self) -> Self::Output {
                Point(self.0.neg())
            }
        }

        impl<'a> std::ops::Neg for &'a Point {
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

        impl<'a, 'b> std::ops::Mul<&'b Scalar> for &'a Point {
            type Output = Point;

            fn mul(self, other: &'b Scalar) -> Point {
                self.scale(other)
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

        impl<'a, 'b> std::ops::Sub<&'b Point> for &'a Point {
            type Output = Point;

            fn sub(self, other: &'b Point) -> Point {
                self + (-other)
            }
        }

        impl std::ops::Sub<Point> for Point {
            type Output = Point;

            fn sub(self, other: Point) -> Point {
                &self - &other
            }
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
