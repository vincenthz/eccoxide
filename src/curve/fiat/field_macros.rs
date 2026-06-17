#[doc(hidden)]
#[macro_export]
macro_rules! fiat_field_common_impl {
    ($(#[$outer:meta])* $FE:ident, $SIZE_BITS:expr, $FE_LIMBS_SIZE:expr, $fiat_constr:ident, $fiat_add:ident, $fiat_sub:ident, $fiat_mul:ident, $fiat_square:ident, $fiat_opp:ident, $fiat_nonzero:ident, $fiat_selectznz:ident) => {
        $(#[$outer])*
        #[derive(Clone)]
        pub struct $FE($fiat_constr);

        impl PartialEq for $FE {
            fn eq(&self, other: &Self) -> bool {
                self.ct_eq(other).is_true()
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

        impl CtZero for $FE {
            fn ct_zero(&self) -> Choice {
                let mut out = 0u64;
                $fiat_nonzero(&mut out, &self.0.0);
                out.ct_zero()
            }
            fn ct_nonzero(&self) -> Choice {
                let mut out = 0u64;
                $fiat_nonzero(&mut out, &self.0.0);
                out.ct_nonzero()
            }
        }

        impl CtEqual<$FE> for $FE {
            fn ct_eq(&self, other: &$FE) -> Choice {
                let r = self - other;
                r.ct_zero()
            }
        }

        impl $FE {
            /// Size in bits of this element of the field
            pub const SIZE_BITS: usize = $SIZE_BITS;

            /// Size in bytes of this element of the field
            pub const SIZE_BYTES: usize = (Self::SIZE_BITS + 7) / 8;

            /// the zero constant (additive identity)
            pub const fn zero() -> Self {
                Self::init([0u64; $FE_LIMBS_SIZE])
            }

            pub fn is_zero(&self) -> bool {
                let mut cond = 0;
                $fiat_nonzero(&mut cond, &self.0.0);
                cond == 0
            }

            /// The one constant (multiplicative identity)
            pub const fn one() -> Self {
                let mut limbs = [0u64; $FE_LIMBS_SIZE];
                limbs[0] = 1;
                Self::init(limbs)
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
            ///
            /// Always true: `self.square() == self * self`
            pub fn square(&self) -> Self {
                let mut out = $fiat_constr([0u64; $FE_LIMBS_SIZE]);
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
            pub const fn double(&self) -> Self {
                let mut out = $fiat_constr([0u64; $FE_LIMBS_SIZE]);
                $fiat_add(&mut out, &self.0, &self.0);
                $FE(out)
            }

            /// Constant-time select: returns `a` if `cond` is true, otherwise `b`.
            ///
            /// Delegates to the fiat-crypto `selectznz` routine, which performs
            /// a branch-free limb-by-limb conditional move. This works whatever
            /// the internal representation is, as a field element and its limbs
            /// are in a one-to-one correspondence.
            pub fn ct_select(cond: Choice, a: &Self, b: &Self) -> Self {
                // selectznz(out, c, x, y) computes `out = if c == 0 { x } else { y }`,
                // so pass `b` then `a` to return `a` when `cond` is true.
                let mut out = $fiat_constr([0u64; $FE_LIMBS_SIZE]);
                $fiat_selectznz(&mut out.0, cond.to_u1(), &b.0 .0, &a.0 .0);
                $FE(out)
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

            /// Similar to 'from_bytes' but take values from a slice.
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

            /// Initialize from a wide buffer of random data.
            ///
            /// The difference with 'from_bytes' or 'from_slice' is that it takes
            /// a random initialized buffer and used modulo operation to initialize
            /// as a field element, but due to inherent bias in modulo operation
            /// we take a double sized buffer.
            ///
            /// The buffer is interpreted as a big-endian integer and reduced
            /// modulo the field characteristic.
            ///
            /// This runs in constant time with respect to the input
            pub fn init_from_wide_bytes(random: [u8; Self::SIZE_BYTES * 2]) -> Self {
                // Reduce the wide big-endian integer modulo p with Horner's
                // method to rewrite the polynomial in nested form: acc = acc * 256 + byte.
                let c256 = Self::from_u64(256);
                let mut acc = Self::zero();
                for b in random.iter() {
                    acc = &acc * &c256 + Self::from_u64(*b as u64);
                }
                acc
            }
        }

        impl std::ops::Neg for $FE {
            type Output = $FE;

            fn neg(self) -> Self::Output {
                let mut out = $fiat_constr([0u64; $FE_LIMBS_SIZE]);
                $fiat_opp(&mut out, &self.0);
                $FE(out)
            }
        }

        impl std::ops::Neg for &$FE {
            type Output = $FE;

            fn neg(self) -> Self::Output {
                let mut out = $fiat_constr([0u64; $FE_LIMBS_SIZE]);
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
                let mut out = $fiat_constr([0u64; $FE_LIMBS_SIZE]);
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
                let mut out = $fiat_constr([0u64; $FE_LIMBS_SIZE]);
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
                let mut out = $fiat_constr([0u64; $FE_LIMBS_SIZE]);
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
            const ZERO : $FE = $FE::zero();
            const ONE : $FE = $FE::one();

            fn is_zero(&self) -> bool {
                self.is_zero()
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
            fn ct_select(cond: Choice, a: &$FE, b: &$FE) -> $FE {
                $FE::ct_select(cond, a, b)
            }
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! fiat_field_montgomery_impl {
    ($(#[$outer:meta])* $FE:ident, $SIZE_BITS:expr, $FIELD_P_LIMBS:expr, $FE_LIMBS_SIZE:expr, $fiat_constr:ident, $fiat_nonzero:ident, $fiat_add:ident, $fiat_sub:ident, $fiat_mul:ident, $fiat_square:ident, $fiat_opp:ident, $fiat_to_bytes:ident, $fiat_from_bytes:ident, $fiat_constr_montgomery:ident, $fiat_to_montgomery:ident, $fiat_from_montgomery:ident, $fiat_selectznz:ident, $fiat_msat:ident, $fiat_divstep:ident, $fiat_divstep_precomp:ident) => {
        crate::fiat_field_common_impl!(
            $(#[$outer])*
            $FE,
            $SIZE_BITS,
            $FE_LIMBS_SIZE,
            $fiat_constr_montgomery,
            $fiat_add,
            $fiat_sub,
            $fiat_mul,
            $fiat_square,
            $fiat_opp,
            $fiat_nonzero,
            $fiat_selectznz
        );

        impl $FE {
            const fn init(current: [u64; $FE_LIMBS_SIZE]) -> Self {
                let mut out = $fiat_constr_montgomery([0u64; $FE_LIMBS_SIZE]);
                $fiat_to_montgomery(&mut out, &$fiat_constr(current));
                Self(out)
            }

            pub const fn from_u64(n: u64) -> Self {
                let mut limbs = [0u64; $FE_LIMBS_SIZE];
                limbs[0] = n;
                Self::init(limbs)
            }


            /// Get the sign of the field element
            pub fn sign(&self) -> Sign {
                let mut out = $fiat_constr([0u64; $FE_LIMBS_SIZE]);
                $fiat_from_montgomery(&mut out, &self.0);
                if out[0] & 1 == 1 {
                    Sign::Negative
                } else {
                    Sign::Positive
                }
            }

            // there's no really negative number in Fp, but if high bit is set ...
            pub fn is_negative(&self) -> bool {
                let mut out = $fiat_constr([0u64; $FE_LIMBS_SIZE]);
                $fiat_from_montgomery(&mut out, &self.0);
                (out[0] & 1) != 0
            }

            /// Initialize a new scalar from its bytes representation (BE) without size check
            ///
            /// This doesn't verify if the element represented fits in the field,
            /// so that run the risk of having elements that are greater or equal
            /// to the order of the field
            pub const fn from_bytes_unchecked(bytes: &[u8; Self::SIZE_BYTES]) -> Self {
                let mut buf = [0u8; Self::SIZE_BYTES];
                buf.copy_from_slice(bytes);
                buf.reverse(); // swap endianness

                let mut out = $fiat_constr([0u64; $FE_LIMBS_SIZE]);
                let mut out_mont = $fiat_constr_montgomery([0u64; $FE_LIMBS_SIZE]);
                $fiat_from_bytes(&mut out.0, &buf);
                $fiat_to_montgomery(&mut out_mont, &out);
                $FE(out_mont)
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

                let mut out = $fiat_constr([0u64; $FE_LIMBS_SIZE]);
                let mut out_mont = $fiat_constr_montgomery([0u64; $FE_LIMBS_SIZE]);
                $fiat_from_bytes(&mut out.0, &buf);

                // modulus limbs, least-significant first
                let mut p = [0u64; $FE_LIMBS_SIZE];
                for (dst, src) in p.iter_mut().zip($FIELD_P_LIMBS.iter().rev()) {
                    *dst = *src;
                }

                // TODO: non constant
                if LimbsLE::ct_lt(LimbsLE(&out.0), LimbsLE(&p[..])).is_true() {
                    $fiat_to_montgomery(&mut out_mont, &out);
                    Some($FE(out_mont))
                } else {
                    None
                }
            }

            /// Output the scalar bytes representation (BE)
            pub const fn to_bytes(&self) -> [u8; Self::SIZE_BYTES] {
                let mut out_normal = $fiat_constr([0u64; $FE_LIMBS_SIZE]);
                let mut out = [0u8; Self::SIZE_BYTES];
                $fiat_from_montgomery(&mut out_normal, &self.0);
                $fiat_to_bytes(&mut out, &out_normal.0);
                out.reverse(); // swap endianness
                out
            }

            /// Get the multiplicative inverse using the Bernstein-Yang
            /// "safegcd" constant-time modular inversion algorithm.
            ///
            /// This is an alternative to the Fermat-little-theorem based
            /// `inverse` (which computes `self^(p-2)` via a curve-specific
            /// addition chain). It builds on the fiat-crypto generated
            /// `msat`/`divstep`/`divstep_precomp` primitives and runs in
            /// constant time with respect to the input.
            ///
            /// Note that fiat-crypto emits a single-step `divstep` (rather
            /// than the batched "jump" variant), so this is typically slower
            /// than `inverse`; it is provided as a representation-agnostic
            /// alternative that does not require a hand-written addition chain.
            ///
            /// Note that 0 doesn't have a multiplicative inverse and will
            /// result in a panic
            pub fn inverse_safegcd(&self) -> Self {
                use crate::mp::ct::CtZero;

                assert!(!self.is_zero());

                // the saturated two's-complement representation of f and g
                // uses one extra limb compared to a field element.
                const SAT_LIMBS: usize = $FE_LIMBS_SIZE + 1;

                // f starts as the modulus m (in saturated form)
                let mut f = [0u64; SAT_LIMBS];
                $fiat_msat(&mut f);

                // The number of divsteps is fixed by the bit length of the
                // modulus, which is public, so this loop bound leaks nothing
                // about the secret input. It must match the exponent baked
                // into the `divstep_precomp` constant.
                let mut len_prime = 0usize;
                let mut i = SAT_LIMBS;
                while i > 0 {
                    i -= 1;
                    if f[i] != 0 {
                        len_prime = i * 64 + (64 - f[i].leading_zeros() as usize);
                        break;
                    }
                }
                let iterations = (49 * len_prime + if len_prime < 46 { 80 } else { 57 }) / 17;

                // g starts as the integer value of self, i.e. taken out of the
                // Montgomery domain and zero-extended into the saturated form.
                let mut a_std = $fiat_constr([0u64; $FE_LIMBS_SIZE]);
                $fiat_from_montgomery(&mut a_std, &self.0);
                let mut g = [0u64; SAT_LIMBS];
                let mut j = 0;
                while j < $FE_LIMBS_SIZE {
                    g[j] = a_std.0[j];
                    j += 1;
                }

                // v and r track the coefficient of `a` for f and g in the
                // Montgomery domain: f = m has coefficient 0 (v = 0) and g = a
                // has coefficient 1 (r = Montgomery one).
                let mut v = [0u64; $FE_LIMBS_SIZE];
                let mut r = (Self::one().0).0;
                let mut d: u64 = 1;

                let mut step = 0;
                while step < iterations {
                    let mut nd = 0u64;
                    let mut nf = [0u64; SAT_LIMBS];
                    let mut ng = [0u64; SAT_LIMBS];
                    let mut nv = [0u64; $FE_LIMBS_SIZE];
                    let mut nr = [0u64; $FE_LIMBS_SIZE];
                    $fiat_divstep(&mut nd, &mut nf, &mut ng, &mut nv, &mut nr, d, &f, &g, &v, &r);
                    d = nd;
                    f = nf;
                    g = ng;
                    v = nv;
                    r = nr;
                    step += 1;
                }

                // After the divsteps, the inverse is in v up to the sign of f:
                // negate v when f ended up negative (top bit of its most
                // significant limb set).
                let v_fe = $FE($fiat_constr_montgomery(v));
                let neg_v = -&v_fe;
                let f_is_negative = (f[SAT_LIMBS - 1] >> 63).ct_nonzero();
                let v_signed = $FE::ct_select(f_is_negative, &neg_v, &v_fe);

                // multiply by the precomputed constant ((m-1)/2)^iterations to
                // undo the factor of 2 accumulated at every divstep.
                let mut precomp = [0u64; $FE_LIMBS_SIZE];
                $fiat_divstep_precomp(&mut precomp);
                let precomp_fe = $FE($fiat_constr_montgomery(precomp));

                &v_signed * &precomp_fe
            }
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! fiat_field_solinas_impl {
    ($(#[$outer:meta])* $FE:ident, $SIZE_BITS:expr, $FIELD_P_BYTES:expr, $FE_LIMBS_SIZE:expr, $fiat_tight_constr:ident, $fiat_loose_constr:ident, $fiat_tighten:ident, $fiat_relax:ident, $fiat_nonzero:ident, $fiat_add:ident, $fiat_sub:ident, $fiat_mul:ident, $fiat_square:ident, $fiat_opp:ident, $fiat_to_bytes:ident, $fiat_from_bytes:ident, $fiat_selectznz:ident) => {
        crate::fiat_field_common_impl!(
            $FE,
            $SIZE_BITS,
            $FE_LIMBS_SIZE,
            $fiat_tight_constr,
            $fiat_add,
            $fiat_sub,
            $fiat_mul,
            $fiat_square,
            $fiat_opp,
            $fiat_nonzero,
            $fiat_selectznz
        );

        impl $FE {
            // build from a little endian limbs size
            //
            // probably should be removed from unsaturated solinas strategy, as this easy
            // to introduce serious bugs..
            const fn init(current: [u64; $FE_LIMBS_SIZE]) -> Self {
                Self($fiat_tight_constr(current))
            }

            pub const fn from_u64(n: u64) -> Self {
                // unsatured solinas run the risk of overflow, so use from_bytes
                // no risk of running into the P limit with a u64
                let [n0, n1, n2, n3, n4, n5, n6, n7] = n.to_be_bytes();
                let mut bytes = [0u8; Self::SIZE_BYTES];
                bytes[Self::SIZE_BYTES - 8] = n0;
                bytes[Self::SIZE_BYTES - 7] = n1;
                bytes[Self::SIZE_BYTES - 6] = n2;
                bytes[Self::SIZE_BYTES - 5] = n3;
                bytes[Self::SIZE_BYTES - 4] = n4;
                bytes[Self::SIZE_BYTES - 3] = n5;
                bytes[Self::SIZE_BYTES - 2] = n6;
                bytes[Self::SIZE_BYTES - 1] = n7;
                //bytes[Self::SIZE_BYTES - 8..].copy_from_slice(&n.to_be_bytes());
                Self::from_bytes_unchecked(&bytes)
            }

            /// Get the sign of the field element
            pub fn sign(&self) -> Sign {
                let out = &self.0;
                if out[0] & 1 == 1 {
                    Sign::Negative
                } else {
                    Sign::Positive
                }
            }

            // there's no really negative number in Fp, but if high bit is set ...
            pub fn is_negative(&self) -> bool {
                let out = &self.0;
                (out[0] & 1) != 0
            }

            /// Initialize a new scalar from its bytes representation (BE)
            ///
            /// This doesn't verify if the element represented fits in the field,
            /// so that run the run of having elements that are greater or equal
            /// than the order of the field
            pub const fn from_bytes_unchecked(bytes: &[u8; Self::SIZE_BYTES]) -> Self {
                let mut buf = [0u8; Self::SIZE_BYTES];
                buf.copy_from_slice(bytes);
                buf.reverse(); // swap endianness

                let mut out = $fiat_tight_constr([0u64; $FE_LIMBS_SIZE]);
                $fiat_from_bytes(&mut out, &buf);
                $FE(out)
            }

            /// Initialize a new scalar from its bytes representation (BE)
            ///
            /// If the represented value overflow the field element size,
            /// then None is returned.
            pub fn from_bytes(bytes: &[u8; Self::SIZE_BYTES]) -> Option<Self> {
                use crate::mp::ct::CtLesser;

                let mut buf = [0u8; Self::SIZE_BYTES];
                buf.copy_from_slice(bytes);
                buf.reverse(); // swap endianness

                let mut out = $fiat_tight_constr([0u64; $FE_LIMBS_SIZE]);
                $fiat_from_bytes(&mut out, &buf);

                // TODO: non constant
                if <&[u8; Self::SIZE_BYTES]>::ct_lt(bytes, &$FIELD_P_BYTES).is_true() {
                    Some($FE(out))
                } else {
                    None
                }
            }

            /// Output the scalar bytes representation (BE)
            pub const fn to_bytes(&self) -> [u8; Self::SIZE_BYTES] {
                let mut out = [0u8; Self::SIZE_BYTES];
                $fiat_to_bytes(&mut out, &self.0);
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
        fn fe_u64(i: u64) -> $FE {
            let mut x8 = [0u8; $FE::SIZE_BYTES];
            let b7 = (i >> 56) as u8;
            let b6 = (i >> 48) as u8;
            let b5 = (i >> 40) as u8;
            let b4 = (i >> 32) as u8;
            let b3 = (i >> 24) as u8;
            let b2 = (i >> 16) as u8;
            let b1 = (i >> 8) as u8;
            let b0 = i as u8;

            x8[$FE::SIZE_BYTES - 8] = b7;
            x8[$FE::SIZE_BYTES - 7] = b6;
            x8[$FE::SIZE_BYTES - 6] = b5;
            x8[$FE::SIZE_BYTES - 5] = b4;
            x8[$FE::SIZE_BYTES - 4] = b3;
            x8[$FE::SIZE_BYTES - 3] = b2;
            x8[$FE::SIZE_BYTES - 2] = b1;
            x8[$FE::SIZE_BYTES - 1] = b0;

            $FE::from_bytes(&x8).expect("working fe_u64")
        }

        fn random_fe_small() -> [$FE; 4] {
            [fe_u64(250), fe_u64(255), fe_u64(256), fe_u64(280)]
        }

        fn eq_small(v1: u64, v2: u64) {
            let f1 = $FE::from_u64(v1);
            let f2 = $FE::from_u64(v2);
            assert_eq!(f1 == f2, v1 == v2, "{} == {}", v1, v2)
        }

        fn ne_small(v1: u64, v2: u64) {
            let f1 = $FE::from_u64(v1);
            let f2 = $FE::from_u64(v2);
            assert_eq!(f1 != f2, v1 != v2, "{} != {}", v1, v2)
        }

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
        fn eq() {
            for i in 1..56 {
                eq_small(i, i);
                eq_small(i, i + 1);
                eq_small(i - 1, i);
                ne_small(i, i);
                ne_small(i, i + 1);
                ne_small(i - 1, i);
            }
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
        fn sub_and_opp_add() {
            for fe in random_fe_small() {
                let zero1 = &fe + (-&fe);
                let zero2 = &fe - &fe;
                assert_eq!(zero1, $FE::zero(), "add opposite not zero");
                assert_eq!(zero2, $FE::zero(), "sub same not zero");
            }
        }

        #[test]
        fn u64_constructor() {
            for i in 250..260 {
                let f1 = $FE::from_u64(i);
                let f2 = fe_u64(i);
                assert_eq!(f1, f2)
            }
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

        #[test]
        fn wide_bytes() {
            assert_eq!(
                $FE::init_from_wide_bytes([0u8; $FE::SIZE_BYTES * 2]),
                $FE::zero()
            );

            let mut wide = [0u8; $FE::SIZE_BYTES * 2];
            wide[$FE::SIZE_BYTES * 2 - 1] = 5;
            assert_eq!($FE::init_from_wide_bytes(wide), $FE::from_u64(5));
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! fiat_field_safegcd_unittest {
    ($FE:ident) => {
        #[test]
        fn inverse_safegcd_matches_fermat() {
            // cross-check the Bernstein-Yang inverse against the (independent)
            // Fermat addition-chain inverse, and the defining property a*a^-1=1.
            for i in 1..200u64 {
                let fe = $FE::from_u64(i);
                let via_flt = fe.inverse();
                let via_by = fe.inverse_safegcd();
                assert_eq!(via_flt, via_by, "safegcd != fermat for {}", i);
                assert_eq!(&fe * &via_by, $FE::one(), "a * a^-1 != 1 for {}", i);
            }

            // also exercise some large/structured values built from bytes
            let mut bytes = [0xa5u8; $FE::SIZE_BYTES];
            bytes[0] = 0; // keep it below the modulus
            let fe = $FE::from_bytes(&bytes).expect("below modulus");
            assert_eq!(fe.inverse(), fe.inverse_safegcd());
            assert_eq!(&fe * &fe.inverse_safegcd(), $FE::one());
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
            let tested = 56;
            for i in 2..tested {
                let f = $FE::from_u64(i);
                match f.sqrt().into_option() {
                    None => println!("{} no sqrt", i),
                    Some(r) => {
                        assert_eq!(&r * &r, f, "$FE returns a sqrt for {} that is not valid", i);
                        found += 1
                    }
                }
            }
            assert!(
                found > 1,
                "not enough sqrt found={} tested={}",
                found,
                tested - 1
            )
        }
    };
}
