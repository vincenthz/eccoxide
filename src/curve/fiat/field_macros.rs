#[doc(hidden)]
#[macro_export]
macro_rules! fiat_field_common_impl {
    ($FE:ident, $SIZE_BITS:expr, $FE_LIMBS_SIZE:expr, $fiat_add:ident, $fiat_sub:ident, $fiat_mul:ident, $fiat_square:ident, $fiat_opp:ident) => {
        #[derive(Clone)]
        pub struct $FE([u64; $FE_LIMBS_SIZE]);

        impl PartialEq for $FE {
            fn eq(&self, other: &Self) -> bool {
                &self.0 == &other.0
            }
        }
        impl Eq for $FE {}

        impl std::fmt::Debug for $FE {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                for b in &self.to_bytes()[..] {
                    write!(f, "{:02x}", b)?
                }
                Ok(())
            }
        }

        impl std::fmt::Display for $FE {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                for b in &self.to_bytes()[..] {
                    write!(f, "{:02x}", b)?
                }
                Ok(())
            }
        }
        impl $FE {
            pub const SIZE_BITS: usize = $SIZE_BITS;
            pub const SIZE_BYTES: usize = (Self::SIZE_BITS + 7) / 8;

            /// the zero constant (additive identity)
            pub fn zero() -> Self {
                Self::init(&[0u64; $FE_LIMBS_SIZE])
            }

            /// The one constant (multiplicative identity)
            pub fn one() -> Self {
                let mut limbs = [0u64; $FE_LIMBS_SIZE];
                limbs[0] = 1;
                Self::init(&limbs)
            }

            pub fn from_u64(n: u64) -> Self {
                let mut limbs = [0u64; $FE_LIMBS_SIZE];
                limbs[0] = n;
                Self::init(&limbs)
            }

            pub fn to_string(&self) -> String {
                let mut s = String::new();
                let bytes = self.to_bytes();
                for b in bytes.iter() {
                    s.push_str(&format!("{:02x}", b));
                }
                s
            }

            /// Return a new element that is the square of this one
            pub fn square(&self) -> Self {
                let mut out = [0u64; $FE_LIMBS_SIZE];
                $fiat_square(&mut out, &self.0);
                Self(out)
            }

            /// Repeadtly square
            fn square_rep(&self, count: usize) -> Self {
                let mut x = self.square();
                for _ in 1..count {
                    x = x.square();
                }
                x
            }

            /// Double the field element, this is equivalent to 2*self or self+self, but can be implemented faster
            pub fn double(&self) -> Self {
                let mut out = [0u64; $FE_LIMBS_SIZE];
                $fiat_add(&mut out, &self.0, &self.0);
                $FE(out)
            }

            pub fn triple(&self) -> Self {
                self + self + self
            }

            /// Compute the field element raised to a power of n, modulus p
            pub fn power_u64(&self, n: u64) -> Self {
                if n == 0 {
                    Self::one()
                } else if n == 1 {
                    self.clone()
                } else if n == 2 {
                    self.square()
                } else {
                    let mut a = self.clone();
                    let mut q = Self::one();

                    for i in 0..64 {
                        if n & (1 << i) != 0 {
                            q = q * &a;
                        }
                        a = a.square();
                    }
                    q
                }
            }

            /// Compute the field element raised to a power of n, modulus p
            pub fn power(&self, limbs: &[u8]) -> Self {
                let mut a = self.clone();
                let mut q = Self::one();

                for limb in limbs.iter().rev() {
                    for i in 0..8 {
                        if limb & (1 << i) != 0 {
                            q = q * &a;
                        }
                        a = a.square();
                    }
                }
                q
            }

            /// Similar to from_bytes but take values from a slice.
            ///
            /// If the slice is not of the right size, then None is returned
            pub fn from_slice(slice: &[u8]) -> Option<Self> {
                if slice.len() != Self::SIZE_BYTES {
                    return None;
                }
                let mut buf = [0u8; Self::SIZE_BYTES];
                buf.copy_from_slice(slice);
                Self::from_bytes(&buf)
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

            // Initialize from a wide buffer of random data.
            //
            // The difference with 'from_bytes' or 'from_slice' is that it takes
            // a random initialized buffer and used modulo operation to initialize
            // as a field element, but due to inherent bias in modulo operation
            // we take a double sized buffer.
            //pub fn init_from_wide_bytes(_random: [u8; Self::SIZE_BYTES * 2]) -> Self {
            //    todo!()
            //}
        }

        impl std::ops::Neg for $FE {
            type Output = $FE;

            fn neg(self) -> Self::Output {
                let mut out = [0u64; $FE_LIMBS_SIZE];
                $fiat_opp(&mut out, &self.0);
                $FE(out)
            }
        }

        impl std::ops::Neg for &$FE {
            type Output = $FE;

            fn neg(self) -> Self::Output {
                let mut out = [0u64; $FE_LIMBS_SIZE];
                $fiat_opp(&mut out, &self.0);
                $FE(out)
            }
        }

        // ****************
        // Scalar Addition
        // ****************

        impl<'a, 'b> std::ops::Add<&'b $FE> for &'a $FE {
            type Output = $FE;

            fn add(self, other: &'b $FE) -> $FE {
                let mut out = [0u64; $FE_LIMBS_SIZE];
                $fiat_add(&mut out, &self.0, &other.0);
                $FE(out)
            }
        }

        impl<'a> std::ops::Add<$FE> for &'a $FE {
            type Output = $FE;

            fn add(self, other: $FE) -> $FE {
                self + &other
            }
        }

        impl<'b> std::ops::Add<&'b $FE> for $FE {
            type Output = $FE;

            fn add(self, other: &'b $FE) -> $FE {
                &self + other
            }
        }

        impl std::ops::Add<$FE> for $FE {
            type Output = $FE;

            fn add(self, other: $FE) -> $FE {
                &self + &other
            }
        }

        // *******************
        // Scalar Subtraction
        // *******************

        impl<'a, 'b> std::ops::Sub<&'b $FE> for &'a $FE {
            type Output = $FE;

            fn sub(self, other: &'b $FE) -> $FE {
                let mut out = [0u64; $FE_LIMBS_SIZE];
                $fiat_sub(&mut out, &self.0, &other.0);
                $FE(out)
            }
        }

        impl<'a> std::ops::Sub<$FE> for &'a $FE {
            type Output = $FE;

            fn sub(self, other: $FE) -> $FE {
                self - &other
            }
        }

        impl<'b> std::ops::Sub<&'b $FE> for $FE {
            type Output = $FE;

            fn sub(self, other: &'b $FE) -> $FE {
                &self - other
            }
        }

        impl std::ops::Sub<$FE> for $FE {
            type Output = $FE;

            fn sub(self, other: $FE) -> $FE {
                &self - &other
            }
        }

        // **********************
        // Scalar Multiplication
        // **********************

        impl<'a, 'b> std::ops::Mul<&'b $FE> for &'a $FE {
            type Output = $FE;

            fn mul(self, other: &'b $FE) -> $FE {
                let mut out = [0u64; $FE_LIMBS_SIZE];
                $fiat_mul(&mut out, &self.0, &other.0);
                $FE(out)
            }
        }

        impl<'b> std::ops::Mul<&'b $FE> for $FE {
            type Output = $FE;

            fn mul(self, other: &'b $FE) -> $FE {
                &self * other
            }
        }

        impl<'a, 'b> std::ops::Mul<$FE> for &'a $FE {
            type Output = $FE;

            fn mul(self, other: $FE) -> $FE {
                self * &other
            }
        }

        impl std::ops::Mul<$FE> for $FE {
            type Output = $FE;

            fn mul(self, other: $FE) -> $FE {
                &self * &other
            }
        }

        impl From<u64> for $FE {
            fn from(v: u64) -> $FE {
                $FE::from_u64(v)
            }
        }

        impl Field for $FE {
            fn zero() -> $FE {
                $FE::zero()
            }
            fn is_zero(&self) -> bool {
                self.is_zero()
            }
            fn one() -> $FE {
                $FE::one()
            }
            fn sign(&self) -> Sign {
                self.sign()
            }
            fn double(&self) -> $FE {
                self.double()
            }
            fn inverse(&self) -> $FE {
                self.inverse()
            }
            fn square(&self) -> $FE {
                self.square()
            }
            fn cube(&self) -> $FE {
                self.square() * self
            }
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! fiat_field_ops_impl {
    ($FE:ident, $SIZE_BITS:expr, $FIELD_P_LIMBS:expr, $FE_LIMBS_SIZE:expr, $fiat_nonzero:ident, $fiat_add:ident, $fiat_sub:ident, $fiat_mul:ident, $fiat_square:ident, $fiat_opp:ident, $fiat_to_bytes:ident, $fiat_from_bytes:ident, $fiat_to_montgomery:ident, $fiat_from_montgomery:ident) => {
        crate::fiat_field_common_impl!(
            $FE,
            $SIZE_BITS,
            $FE_LIMBS_SIZE,
            $fiat_add,
            $fiat_sub,
            $fiat_mul,
            $fiat_square,
            $fiat_opp
        );

        impl CtZero for &$FE {
            fn ct_zero(f: &$FE) -> Choice {
                let mut out = 0;
                $fiat_nonzero(&mut out, &f.0);
                CtZero::ct_zero(out)
            }
            fn ct_nonzero(f: &$FE) -> Choice {
                let mut out = 0;
                $fiat_nonzero(&mut out, &f.0);
                CtZero::ct_nonzero(out)
            }
        }

        impl CtEqual<$FE> for $FE {
            fn ct_eq(&self, other: &$FE) -> Choice {
                let r = self - other;
                CtZero::ct_zero(&r)
            }
        }

        impl $FE {
            fn init(current: &[u64; $FE_LIMBS_SIZE]) -> Self {
                let mut out = [0u64; $FE_LIMBS_SIZE];
                $fiat_to_montgomery(&mut out, &current);
                Self(out)
            }

            pub fn is_zero(&self) -> bool {
                let mut cond = 0;
                $fiat_nonzero(&mut cond, &self.0);
                cond == 0
            }

            /// Get the sign of the field element
            pub fn sign(&self) -> Sign {
                let mut out = [0u64; $FE_LIMBS_SIZE];
                $fiat_from_montgomery(&mut out, &self.0);
                if out[0] & 1 == 1 {
                    Sign::Negative
                } else {
                    Sign::Positive
                }
            }

            // there's no really negative number in Fp, but if high bit is set ...
            pub fn is_negative(&self) -> bool {
                let mut out = [0u64; $FE_LIMBS_SIZE];
                $fiat_from_montgomery(&mut out, &self.0);
                (out[0] & 1) != 0
            }

            /// Initialize a new scalar from its bytes representation (BE)
            ///
            /// If the represented value overflow the field element size,
            /// then None is returned.
            pub fn from_bytes(bytes: &[u8; Self::SIZE_BYTES]) -> Option<Self> {
                use crate::mp::ct::CtLesser;
                use crate::mp::limbs::LimbsLE;

                let mut buf = [0u8; Self::SIZE_BYTES];
                buf.copy_from_slice(bytes);
                buf.reverse(); // swap endianness

                let mut out = [0u64; $FE_LIMBS_SIZE];
                let mut out_mont = [0u64; $FE_LIMBS_SIZE];
                $fiat_from_bytes(&mut out, &buf);

                let p = $FIELD_P_LIMBS.iter().rev().copied().collect::<Vec<_>>();

                // TODO: non constant
                if LimbsLE::ct_lt(LimbsLE(&out), LimbsLE(&p[..])).is_true() {
                    $fiat_to_montgomery(&mut out_mont, &out);
                    Some($FE(out_mont))
                } else {
                    None
                }
            }

            /// Output the scalar bytes representation (BE)
            pub fn to_bytes(&self) -> [u8; Self::SIZE_BYTES] {
                let mut out_normal = [0u64; $FE_LIMBS_SIZE];
                let mut out = [0u8; Self::SIZE_BYTES];
                $fiat_from_montgomery(&mut out_normal, &self.0);
                $fiat_to_bytes(&mut out, &out_normal);
                out.reverse(); // swap endianness
                out
            }
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! fiat_field_sqrt_define {
    ($FE:ident) => {
        impl FieldSqrt for $FE {
            fn sqrt(&self) -> CtOption<$FE> {
                self.sqrt()
            }
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! fiat_field_unittest {
    ($FE:ident) => {
        fn add_small(v1: u64, v2: u64) {
            let f1 = $FE::from_u64(v1);
            let f2 = $FE::from_u64(v2);
            let fr = $FE::from_u64(v1 + v2);
            assert_eq!(f1 + f2, fr)
        }

        fn square_small(v1: u64) {
            let f1 = $FE::from_u64(v1);
            let fr = $FE::from_u64(v1 * v1);
            assert_eq!(f1.square(), fr)
        }

        fn mul_small(v1: u64, v2: u64) {
            let f1 = $FE::from_u64(v1);
            let f2 = $FE::from_u64(v2);
            let fr = $FE::from_u64(v1 * v2);
            assert_eq!(f1 * f2, fr)
        }

        fn power_small(v1: u64, v2: u32) {
            let f1 = $FE::from_u64(v1);
            let fr = $FE::from_u64(v1.pow(v2));
            assert_eq!(f1.power_u64(v2 as u64), fr)
        }

        #[test]
        fn add() {
            add_small(3, 24);
            add_small(0xff01, 1);
            add_small(0x10001, 0x100);
        }

        #[test]
        fn square() {
            square_small(3);
            square_small(0xff01);
            square_small(0x10001);
        }

        #[test]
        fn mul() {
            mul_small(3, 24);
            mul_small(0x0, 1);
            mul_small(0xff01, 1);
            mul_small(0x10001, 0x100);
        }

        #[test]
        fn power() {
            power_small(3, 24);
            power_small(0x2, 9);
            power_small(0xff01, 4);
            power_small(0x13, 13);
        }

        #[test]
        fn sub() {
            let f1 = $FE::from_u64(49);
            let f2 = $FE::from_u64(24);
            let fr = $FE::from_u64(25);

            assert_eq!(f1 - f2, fr)
        }

        #[test]
        fn inverse() {
            for i in 1..124 {
                let fe = $FE::from_u64(i);
                let r = &fe * fe.inverse();
                println!("{} * {} = {}", fe, fe.inverse(), r);
                assert_eq!($FE::one(), r);
            }
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! fiat_field_sqrt_unittest {
    ($FE:ident) => {
        #[test]
        fn sqrt() {
            let mut found = 0;
            for i in 2..34 {
                let f = $FE::from_u64(i);
                match f.sqrt().into_option() {
                    None => println!("{} no sqrt", i),
                    Some(r) => {
                        assert_eq!(&r * &r, f, "$FE returns a sqrt for {} that is not valid", i);
                        found += 1
                    }
                }
            }
            assert!(found > 1, "not enough sqrt found")
        }
    };
}
