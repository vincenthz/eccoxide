use crate::curve::fiat::secp256k1_scalar_64::*;
//use crate::params::sec2::p256k1::PM2_BYTES;

const GM_LIMBS_SIZE: usize = 4;

fn swap_endian(buf: &mut [u8; 32]) {
    for i in 0..16 {
        let v = buf[i];
        buf[i] = buf[31 - i];
        buf[31 - i] = v
    }
}

#[derive(Clone)]
pub struct Scalar([u64; GM_LIMBS_SIZE]);

impl PartialEq for Scalar {
    fn eq(&self, other: &Self) -> bool {
        &self.0 == &other.0
    }
}

impl std::fmt::Debug for Scalar {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for b in &self.to_bytes() {
            write!(f, "{:02x}", b)?
        }
        Ok(())
    }
}

impl Eq for Scalar {}

impl Scalar {
    pub const SIZE_BITS: usize = 256;
    pub const SIZE_BYTES: usize = (Self::SIZE_BITS + 7) / 8;

    /*
    /// init from limbs to internal representation (montgomery)
    fn init(current: &[u64; GM_LIMBS_SIZE]) -> Self {
        let mut out = [0u64; GM_LIMBS_SIZE];
        let mut current_swapped = [0u64; GM_LIMBS_SIZE];
        current_swapped[0] = current[3];
        current_swapped[1] = current[2];
        current_swapped[2] = current[1];
        current_swapped[3] = current[0];
        fiat_secp256k1_scalar_to_montgomery(&mut out, &current_swapped);
        Self(out)
    }
    */

    /// the zero constant (additive identity)
    pub fn zero() -> Self {
        Self::from_bytes(&[
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0,
        ])
        .unwrap()
    }

    /// The one constant (multiplicative identity)
    pub fn one() -> Self {
        Self::from_bytes(&[
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1,
        ])
        .unwrap()
    }

    pub fn from_u64(n: u64) -> Self {
        let mut bytes = [0u8; 32];
        bytes[24..32].copy_from_slice(&n.to_be_bytes());
        Self::from_bytes(&bytes).unwrap()
    }

    pub fn is_zero(&self) -> bool {
        let mut cond = 0;
        fiat_secp256k1_scalar_nonzero(&mut cond, &self.0);
        cond == 0
    }

    // there's no really negative number in Fp, but if high bit is set ...
    pub fn high_bit_set(&self) -> bool {
        todo!()
    }

    /// Self add another Scalar
    pub fn add_assign(&mut self, _other: &Self) {
        todo!()
    }

    pub fn to_string(&self) -> String {
        let mut s = String::new();
        let bytes = self.to_bytes();
        for b in bytes.iter() {
            s.push_str(&format!("{:02x}", b));
        }
        s
    }

    pub fn square(&self) -> Self {
        let mut out = [0u64; GM_LIMBS_SIZE];
        fiat_secp256k1_scalar_square(&mut out, &self.0);
        Self(out)
    }

    fn square_rep(&self, count: usize) -> Self {
        let mut x = self.square();
        for _ in 1..count {
            x = x.square();
        }
        x
    }

    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse
    pub fn inverse(&self) -> Self {
        let x = self;
        let u2 = self.square();
        let x2 = &u2 * x;
        let u5 = &u2 * &x2;
        let x3 = &u5 * &u2;
        let u9 = &x3 * &u2;
        let u11 = &u9 * &u2;
        let u13 = &u11 * &u2;

        let x6 = u13.square().square() * &u11;
        let x8 = x6.square().square() * &x2;

        let x14 = x8.square_rep(6) * &x6;
        let x28 = x14.square_rep(14) * &x14;
        let x56 = x28.square_rep(28) * &x28;
        let x112 = x56.square_rep(56) * &x56;
        let x126 = x112.square_rep(14) * &x14;

        let mut t = x126.square_rep(3) * &u5;
        t = t.square_rep(4) * &x3;
        t = t.square_rep(4) * &u5;
        t = t.square_rep(5) * &u11;
        t = t.square_rep(4) * &u11;
        t = t.square_rep(4) * &x3;
        t = t.square_rep(5) * &x3;
        t = t.square_rep(6) * &u13;
        t = t.square_rep(4) * &u5;
        t = t.square_rep(3) * &x3;
        t = t.square_rep(5) * &u9;

        t = t.square_rep(6) * &u5;
        t = t.square_rep(10) * &x3;
        t = t.square_rep(4) * &x3;
        t = t.square_rep(9) * &x8;
        t = t.square_rep(5) * &u9;
        t = t.square_rep(6) * &u11;
        t = t.square_rep(4) * &u13;
        t = t.square_rep(5) * &x2;
        t = t.square_rep(6) * &u13;
        t = t.square_rep(10) * &u13;
        t = t.square_rep(4) * &u9;
        t = t.square_rep(6) * x;
        t = t.square_rep(8) * &x6;

        t
    }

    /// Double the field element, this is equivalent to 2*self or self+self, but can be implemented faster
    pub fn double(&self) -> Self {
        let mut out = [0u64; GM_LIMBS_SIZE];
        fiat_secp256k1_scalar_add(&mut out, &self.0, &self.0);
        Scalar(out)
    }

    pub fn triple(&self) -> Self {
        let mut out = [0u64; GM_LIMBS_SIZE];
        let mut out2 = [0u64; GM_LIMBS_SIZE];
        fiat_secp256k1_scalar_add(&mut out, &self.0, &self.0);
        fiat_secp256k1_scalar_add(&mut out2, &out, &self.0);
        Scalar(out2)
    }

    pub fn power_(&self, limbs: &[u8]) -> Self {
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

    /// Compute the field element raised to a power of n, modulus p
    pub fn power(&self, n: u64) -> Self {
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

    /// Initialize a new scalar from its bytes representation
    ///
    /// If the represented value overflow the field element size,
    /// then None is returned.
    pub fn from_bytes(bytes: &[u8; Self::SIZE_BYTES]) -> Option<Self> {
        let mut buf = [0u8; Self::SIZE_BYTES];
        buf.copy_from_slice(bytes);
        swap_endian(&mut buf);

        let mut out = [0u64; GM_LIMBS_SIZE];
        let mut out_mont = [0u64; GM_LIMBS_SIZE];
        fiat_secp256k1_scalar_from_bytes(&mut out, &buf);
        fiat_secp256k1_scalar_to_montgomery(&mut out_mont, &out);
        Some(Scalar(out_mont))
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
        swap_endian(&mut buf);

        let mut out = [0u64; GM_LIMBS_SIZE];
        fiat_secp256k1_scalar_from_bytes(&mut out, &buf);
        Some(Scalar(out))
    }

    /// Output the scalar bytes representation
    pub fn to_bytes(&self) -> [u8; Self::SIZE_BYTES] {
        let mut out_normal = [0u64; GM_LIMBS_SIZE];
        let mut out = [0u8; Self::SIZE_BYTES];
        fiat_secp256k1_scalar_from_montgomery(&mut out_normal, &self.0);
        fiat_secp256k1_scalar_to_bytes(&mut out, &out_normal);
        swap_endian(&mut out);
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
    pub fn init_from_wide_bytes(_random: [u8; Self::SIZE_BYTES * 2]) -> Self {
        todo!()
    }
}

impl std::ops::Neg for Scalar {
    type Output = Scalar;

    fn neg(self) -> Self::Output {
        let mut out = [0u64; GM_LIMBS_SIZE];
        fiat_secp256k1_scalar_opp(&mut out, &self.0);
        Scalar(out)
    }
}

impl std::ops::Neg for &Scalar {
    type Output = Scalar;

    fn neg(self) -> Self::Output {
        let mut out = [0u64; GM_LIMBS_SIZE];
        fiat_secp256k1_scalar_opp(&mut out, &self.0);
        Scalar(out)
    }
}

// ****************
// Scalar Addition
// ****************

impl<'a, 'b> std::ops::Add<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    fn add(self, other: &'b Scalar) -> Scalar {
        let mut out = [0u64; GM_LIMBS_SIZE];
        fiat_secp256k1_scalar_add(&mut out, &self.0, &other.0);
        Scalar(out)
    }
}

impl<'a> std::ops::Add<Scalar> for &'a Scalar {
    type Output = Scalar;

    fn add(self, other: Scalar) -> Scalar {
        self + &other
    }
}

impl<'b> std::ops::Add<&'b Scalar> for Scalar {
    type Output = Scalar;

    fn add(self, other: &'b Scalar) -> Scalar {
        &self + other
    }
}

impl std::ops::Add<Scalar> for Scalar {
    type Output = Scalar;

    fn add(self, other: Scalar) -> Scalar {
        &self + &other
    }
}

// *******************
// Scalar Subtraction
// *******************

impl<'a, 'b> std::ops::Sub<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    fn sub(self, other: &'b Scalar) -> Scalar {
        let mut out = [0u64; GM_LIMBS_SIZE];
        fiat_secp256k1_scalar_sub(&mut out, &self.0, &other.0);
        Scalar(out)
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
        let mut out = [0u64; GM_LIMBS_SIZE];
        fiat_secp256k1_scalar_mul(&mut out, &self.0, &other.0);
        Scalar(out)
    }
}

impl<'b> std::ops::Mul<&'b Scalar> for Scalar {
    type Output = Scalar;

    fn mul(self, other: &'b Scalar) -> Scalar {
        &self * other
    }
}

impl<'a, 'b> std::ops::Mul<Scalar> for &'a Scalar {
    type Output = Scalar;

    fn mul(self, other: Scalar) -> Scalar {
        self * &other
    }
}

impl std::ops::Mul<Scalar> for Scalar {
    type Output = Scalar;

    fn mul(self, other: Scalar) -> Scalar {
        &self * &other
    }
}

#[cfg(test)]
mod tests {
    use super::Scalar;

    fn add_small(v1: u64, v2: u64) {
        let f1 = Scalar::from_u64(v1);
        let f2 = Scalar::from_u64(v2);
        let fr = Scalar::from_u64(v1 + v2);
        assert_eq!(f1 + f2, fr)
    }

    fn square_small(v1: u64) {
        let f1 = Scalar::from_u64(v1);
        let fr = Scalar::from_u64(v1 * v1);
        assert_eq!(f1.square(), fr)
    }

    fn mul_small(v1: u64, v2: u64) {
        let f1 = Scalar::from_u64(v1);
        let f2 = Scalar::from_u64(v2);
        let fr = Scalar::from_u64(v1 * v2);
        assert_eq!(f1 * f2, fr)
    }

    fn power_small(v1: u64, v2: u32) {
        let f1 = Scalar::from_u64(v1);
        let fr = Scalar::from_u64(v1.pow(v2));
        assert_eq!(f1.power(v2 as u64), fr)
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
        let f1 = Scalar::from_u64(49);
        let f2 = Scalar::from_u64(24);
        let fr = Scalar::from_u64(25);

        assert_eq!(f1 - f2, fr)
    }

    #[test]
    fn inverse() {
        for i in 1..124 {
            println!("{}", i);
            let fe = Scalar::from_u64(i);
            let r = &fe * fe.inverse();
            assert_eq!(Scalar::one(), r);
        }
    }
}
