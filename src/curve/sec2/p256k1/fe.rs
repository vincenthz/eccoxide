use crate::curve::fiat::secp256k1_64::*;
use crate::mp::ct::{Choice, CtEqual, CtOption, CtZero};
use crate::params::sec2::p256k1::P_LIMBS;

const FE_LIMBS_SIZE: usize = 4;

fn swap_endian(buf: &mut [u8; 32]) {
    for i in 0..16 {
        let v = buf[i];
        buf[i] = buf[31 - i];
        buf[31 - i] = v
    }
}

#[derive(Clone)]
pub struct FieldElement([u64; FE_LIMBS_SIZE]);

impl PartialEq for FieldElement {
    fn eq(&self, other: &Self) -> bool {
        &self.0 == &other.0
    }
}

impl std::fmt::Debug for FieldElement {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for b in &self.to_bytes() {
            write!(f, "{:02x}", b)?
        }
        Ok(())
    }
}

impl CtZero for &FieldElement {
    fn ct_zero(f: &FieldElement) -> Choice {
        let mut out = 0;
        fiat_secp256k1_nonzero(&mut out, &f.0);
        CtZero::ct_zero(out)
    }
    fn ct_nonzero(f: &FieldElement) -> Choice {
        let mut out = 0;
        fiat_secp256k1_nonzero(&mut out, &f.0);
        CtZero::ct_nonzero(out)
    }
}
impl CtEqual for &FieldElement {
    fn ct_eq(f1: &FieldElement, f2: &FieldElement) -> Choice {
        let r = f1 - f2;
        CtZero::ct_zero(&r)
    }
}
impl Eq for FieldElement {}

impl FieldElement {
    pub const SIZE_BITS: usize = 256;
    pub const SIZE_BYTES: usize = (Self::SIZE_BITS + 7) / 8;

    /*
    /// init from limbs to internal representation (montgomery)
    pub(super) fn init(current: &[u64; FE_LIMBS_SIZE]) -> Self {
        let mut out = [0u64; FE_LIMBS_SIZE];
        let mut current_swapped = [0u64; FE_LIMBS_SIZE];
        current_swapped[0] = u64::from_be(current[3]);
        current_swapped[1] = u64::from_be(current[2]);
        current_swapped[2] = u64::from_be(current[1]);
        current_swapped[3] = u64::from_be(current[0]);
        fiat_secp256k1_to_montgomery(&mut out, &current_swapped);
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
        //Self::init(&[0, 0, 0, 0])
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
        fiat_secp256k1_nonzero(&mut cond, &self.0);
        cond == 0
    }

    // there's no really negative number in Fp, but if high bit is set ...
    pub fn high_bit_set(&self) -> bool {
        let mut out = [0u64; FE_LIMBS_SIZE];
        fiat_secp256k1_from_montgomery(&mut out, &self.0);
        (out[0] & 1) != 0
    }

    /// Self add another Scalar
    pub fn add_assign(&mut self, _other: &Self) {
        todo!()
    }

    pub fn square(&self) -> Self {
        let mut out = [0u64; FE_LIMBS_SIZE];
        fiat_secp256k1_square(&mut out, &self.0);
        Self(out)
    }

    fn square_rep(&self, count: usize) -> Self {
        let mut x = self.square();
        for _ in 1..count {
            x = x.square();
        }
        x
    }

    pub fn to_string(&self) -> String {
        let mut s = String::new();
        let bytes = self.to_bytes();
        for b in bytes.iter() {
            s.push_str(&format!("{:02x}", b));
        }
        s
    }

    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse
    fn r_inverse(&self) -> Self {
        let x2 = self.square() * self;
        let x3 = x2.square() * self;
        let x6 = x3.square_rep(3) * &x3;
        let x9 = x6.square_rep(3) * &x3;
        let x11 = x9.square_rep(2) * &x2;
        let x22 = x11.square_rep(11) * &x11;
        let x44 = x22.square_rep(22) * &x22;
        let x88 = x44.square_rep(44) * &x44;
        let x176 = x88.square_rep(88) * &x88;
        let x220 = x176.square_rep(44) * &x44;
        let x223 = x220.square_rep(3) * &x3;

        let mut t1 = x223.square_rep(23) * &x22;
        t1 = t1.square_rep(5) * self;
        t1 = t1.square_rep(3) * &x2;
        t1 = t1.square_rep(2);

        t1 * self
    }

    pub fn inverse(&self) -> Self {
        assert!(!self.is_zero());
        self.r_inverse()
    }

    /// Double the field element, this is equivalent to 2*self or self+self, but can be implemented faster
    pub fn double(&self) -> Self {
        self + self
    }

    pub fn triple(&self) -> Self {
        self + self + self
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
            let mut q = Self::zero();

            for i in 0..64 {
                if n & (1 << i) != 0 {
                    q = q * &a;
                }
                a = a.square();
            }
            q
        }
    }

    /// Compute the square root 'x' of the field element such that x*x = self
    pub fn sqrt(&self) -> CtOption<Self> {
        let x2 = self.square() * self;
        let x3 = x2.square() * self;
        let x6 = x3.square_rep(3) * &x3;
        let x9 = x6.square_rep(3) * &x3;
        let x11 = x9.square_rep(2) * &x2;
        let x22 = x11.square_rep(11) * &x11;
        let x44 = x22.square_rep(22) * &x22;
        let x88 = x44.square_rep(44) * &x44;
        let x176 = x88.square_rep(88) * &x88;
        let x220 = x176.square_rep(44) * &x44;
        let x223 = x220.square_rep(3) * &x3;

        let mut t1 = x223.square_rep(23) * &x22;
        t1 = t1.square_rep(6) * &x2;
        t1 = &t1 * &t1;

        let r = &t1 * &t1;
        let r2 = &r * &r;
        CtOption::from((CtEqual::ct_eq(&r2, self), r))
    }

    /// Initialize a new scalar from its bytes representation
    ///
    /// If the represented value overflow the field element size,
    /// then None is returned.
    pub fn from_bytes(bytes: &[u8; Self::SIZE_BYTES]) -> Option<Self> {
        use crate::mp::ct::CtLesser;
        use crate::mp::limbs::LimbsLE;

        let mut buf = [0u8; Self::SIZE_BYTES];
        buf.copy_from_slice(bytes);
        swap_endian(&mut buf);

        let mut out = [0u64; FE_LIMBS_SIZE];
        let mut out_mont = [0u64; FE_LIMBS_SIZE];
        fiat_secp256k1_from_bytes(&mut out, &buf);

        let p = P_LIMBS.iter().rev().copied().collect::<Vec<_>>();

        if LimbsLE::ct_lt(LimbsLE(&out), LimbsLE(&p[..])).is_true() {
            fiat_secp256k1_to_montgomery(&mut out_mont, &out);
            Some(FieldElement(out_mont))
        } else {
            None
        }
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

    /// Output the scalar bytes representation
    pub fn to_bytes(&self) -> [u8; Self::SIZE_BYTES] {
        let mut out_normal = [0u64; FE_LIMBS_SIZE];
        let mut out = [0u8; Self::SIZE_BYTES];
        fiat_secp256k1_from_montgomery(&mut out_normal, &self.0);
        fiat_secp256k1_to_bytes(&mut out, &out_normal);
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

impl std::ops::Neg for FieldElement {
    type Output = FieldElement;

    fn neg(self) -> Self::Output {
        let mut out = [0u64; FE_LIMBS_SIZE];
        fiat_secp256k1_opp(&mut out, &self.0);
        FieldElement(out)
    }
}

impl std::ops::Neg for &FieldElement {
    type Output = FieldElement;

    fn neg(self) -> Self::Output {
        let mut out = [0u64; FE_LIMBS_SIZE];
        fiat_secp256k1_opp(&mut out, &self.0);
        FieldElement(out)
    }
}

// ****************
// Scalar Addition
// ****************

impl<'a, 'b> std::ops::Add<&'b FieldElement> for &'a FieldElement {
    type Output = FieldElement;

    fn add(self, other: &'b FieldElement) -> FieldElement {
        let mut out = [0u64; FE_LIMBS_SIZE];
        fiat_secp256k1_add(&mut out, &self.0, &other.0);
        FieldElement(out)
    }
}

impl<'a> std::ops::Add<FieldElement> for &'a FieldElement {
    type Output = FieldElement;

    fn add(self, other: FieldElement) -> FieldElement {
        self + &other
    }
}

impl<'b> std::ops::Add<&'b FieldElement> for FieldElement {
    type Output = FieldElement;

    fn add(self, other: &'b FieldElement) -> FieldElement {
        &self + other
    }
}

impl std::ops::Add<FieldElement> for FieldElement {
    type Output = FieldElement;

    fn add(self, other: FieldElement) -> FieldElement {
        &self + &other
    }
}

// *******************
// Scalar Subtraction
// *******************

impl<'a, 'b> std::ops::Sub<&'b FieldElement> for &'a FieldElement {
    type Output = FieldElement;

    fn sub(self, other: &'b FieldElement) -> FieldElement {
        let mut out = [0u64; FE_LIMBS_SIZE];
        fiat_secp256k1_sub(&mut out, &self.0, &other.0);
        FieldElement(out)
    }
}

impl<'a> std::ops::Sub<FieldElement> for &'a FieldElement {
    type Output = FieldElement;

    fn sub(self, other: FieldElement) -> FieldElement {
        self - &other
    }
}

impl<'b> std::ops::Sub<&'b FieldElement> for FieldElement {
    type Output = FieldElement;

    fn sub(self, other: &'b FieldElement) -> FieldElement {
        &self - other
    }
}

impl std::ops::Sub<FieldElement> for FieldElement {
    type Output = FieldElement;

    fn sub(self, other: FieldElement) -> FieldElement {
        &self - &other
    }
}

// **********************
// Scalar Multiplication
// **********************

impl<'a, 'b> std::ops::Mul<&'b FieldElement> for &'a FieldElement {
    type Output = FieldElement;

    fn mul(self, other: &'b FieldElement) -> FieldElement {
        let mut out = [0u64; FE_LIMBS_SIZE];
        fiat_secp256k1_mul(&mut out, &self.0, &other.0);
        FieldElement(out)
    }
}

impl<'b> std::ops::Mul<&'b FieldElement> for FieldElement {
    type Output = FieldElement;

    fn mul(self, other: &'b FieldElement) -> FieldElement {
        &self * other
    }
}

impl<'a, 'b> std::ops::Mul<FieldElement> for &'a FieldElement {
    type Output = FieldElement;

    fn mul(self, other: FieldElement) -> FieldElement {
        self * &other
    }
}

impl std::ops::Mul<FieldElement> for FieldElement {
    type Output = FieldElement;

    fn mul(self, other: FieldElement) -> FieldElement {
        &self * &other
    }
}

#[cfg(test)]
mod tests {
    use super::FieldElement;

    fn add_small(v1: u64, v2: u64) {
        let f1 = FieldElement::from_u64(v1);
        let f2 = FieldElement::from_u64(v2);
        let fr = FieldElement::from_u64(v1 + v2);
        assert_eq!(f1 + f2, fr)
    }

    fn square_small(v1: u64) {
        let f1 = FieldElement::from_u64(v1);
        let fr = FieldElement::from_u64(v1 * v1);
        assert_eq!(f1.square(), fr)
    }

    fn mul_small(v1: u64, v2: u64) {
        let f1 = FieldElement::from_u64(v1);
        let f2 = FieldElement::from_u64(v2);
        let fr = FieldElement::from_u64(v1 * v2);
        assert_eq!(f1 * f2, fr)
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
    fn sub() {
        let f1 = FieldElement::from_u64(49);
        let f2 = FieldElement::from_u64(24);
        let fr = FieldElement::from_u64(25);

        assert_eq!(f1 - f2, fr)
    }

    #[test]
    fn sqrt() {
        for i in 2..34 {
            let f = FieldElement::from_u64(i);
            match f.sqrt().into_option() {
                None => (),
                Some(r) => assert_eq!(
                    &r * &r,
                    f,
                    "FieldElement returns a sqrt for {} that is not valid",
                    i
                ),
            }
        }
    }

    #[test]
    fn inverse() {
        for i in 2..124 {
            println!("{}", i);
            let fe = FieldElement::from_u64(i);
            let r = &fe * fe.inverse();
            assert_eq!(FieldElement::one(), r);
        }
    }
}
