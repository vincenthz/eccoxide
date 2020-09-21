use crate::curve::fiat::secp256k1_scalar_64::*;
use crate::curve::field::{Field, Sign};
use crate::fiat_field_ops_impl;
use crate::mp::ct::{Choice, CtEqual, CtZero};
use crate::params::sec2::p256k1::ORDER_LIMBS;

const GM_LIMBS_SIZE: usize = 4;

fiat_field_ops_impl!(
    Scalar,
    256,
    ORDER_LIMBS,
    GM_LIMBS_SIZE,
    fiat_secp256k1_scalar_nonzero,
    fiat_secp256k1_scalar_to_montgomery,
    fiat_secp256k1_scalar_from_montgomery,
    fiat_secp256k1_scalar_add,
    fiat_secp256k1_scalar_sub,
    fiat_secp256k1_scalar_mul,
    fiat_secp256k1_scalar_square,
    fiat_secp256k1_scalar_opp,
    fiat_secp256k1_scalar_to_bytes,
    fiat_secp256k1_scalar_from_bytes
);

impl Scalar {
    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    /// TODO this will change to being a method of NonZeroScalar
    pub fn inverse(&self) -> Self {
        assert!(!self.is_zero());
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
