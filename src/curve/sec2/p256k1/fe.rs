use crate::curve::fiat::secp256k1_64::*;
use crate::curve::field::{Field, FieldSqrt, Sign};
use crate::mp::ct::{Choice, CtEqual, CtOption, CtZero};
use crate::params::sec2::p256k1::P_LIMBS;
use crate::{fiat_field_ops_impl, fiat_field_sqrt_define};

const FE_LIMBS_SIZE: usize = 4;

fiat_field_ops_impl!(
    FieldElement,
    256,
    P_LIMBS,
    FE_LIMBS_SIZE,
    fiat_secp256k1_nonzero,
    fiat_secp256k1_to_montgomery,
    fiat_secp256k1_from_montgomery,
    fiat_secp256k1_add,
    fiat_secp256k1_sub,
    fiat_secp256k1_mul,
    fiat_secp256k1_square,
    fiat_secp256k1_opp,
    fiat_secp256k1_to_bytes,
    fiat_secp256k1_from_bytes
);
fiat_field_sqrt_define!(FieldElement);

impl FieldElement {
    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    pub fn inverse(&self) -> Self {
        assert!(!self.is_zero());
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

    fn power_small(v1: u64, v2: u32) {
        let f1 = FieldElement::from_u64(v1);
        let fr = FieldElement::from_u64(v1.pow(v2));
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
        for i in 1..124 {
            println!("{}", i);
            let fe = FieldElement::from_u64(i);
            let r = &fe * fe.inverse();
            assert_eq!(FieldElement::one(), r);
        }
    }
}
