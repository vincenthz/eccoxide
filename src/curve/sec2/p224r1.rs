//! Curve p224r1 as defined over the prime field of order 2^224 - 2^96 + 1
use crate::curve::fiat::p224r1_64::*;
use crate::curve::fiat::p224r1_scalar_64::*;
use crate::curve::field::{Field, FieldSqrt, Sign};
use crate::curve::{affine, projective, weierstrass::WeierstrassCurve};
use crate::mp::ct::{Choice, CtEqual, CtOption, CtZero};
use crate::params::sec2::p224r1::*;
use crate::{fiat_define_weierstrass_curve, fiat_define_weierstrass_points};
use crate::{fiat_field_ops_impl, fiat_field_sqrt_define};

const GM_LIMBS_SIZE: usize = 4;
const FE_LIMBS_SIZE: usize = 4;

fiat_field_ops_impl!(
    #[doc = "Element of the prime field Fp where p = 2^224 - 2^96 + 1"]
    FieldElement,
    224,
    P_LIMBS,
    FE_LIMBS_SIZE,
    fiat_p224r1_nonzero,
    fiat_p224r1_add,
    fiat_p224r1_sub,
    fiat_p224r1_mul,
    fiat_p224r1_square,
    fiat_p224r1_opp,
    fiat_p224r1_to_bytes,
    fiat_p224r1_from_bytes,
    montgomery {
        fiat_p224r1_to_montgomery,
        fiat_p224r1_from_montgomery
    }
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
        let x4 = x3.square() * self;
        let x6 = x4.square_rep(2) * &x2;
        let x7 = x6.square() * self;
        let x14 = x7.square_rep(7) * &x7;
        let x28 = x14.square_rep(14) * &x14;
        let x56 = x28.square_rep(28) * &x28;
        let x62 = x56.square_rep(6) * &x6;
        let x90 = x62.square_rep(28) * &x28;
        let x96 = x90.square_rep(6) * &x6;
        let x124 = x96.square_rep(28) * &x28;
        let x127 = x124.square_rep(3) * &x3;
        // 1*127,0*1,1*96
        let t1 = x127.square_rep(97) * &x96; // 1*127 0*1 1*96
        t1
    }

    pub fn is_quadratic_residue(&self) -> Choice {
        let euler = {
            let x2 = self.square() * self;
            let x4 = x2.square_rep(2) * &x2;
            let x8 = x4.square_rep(4) * &x4;
            let x16 = x8.square_rep(8) * &x8;
            let x32 = x16.square_rep(16) * &x16;
            let x64 = x32.square_rep(32) * &x32;
            let x128 = x64.square_rep(64) * &x64;
            // [1*128,0*95]
            x128.square_rep(95)
        };
        euler.ct_eq(&FieldElement::one())
    }

    fn square_add_rep(&self, count: usize) -> Self {
        let mut a = self.clone();
        let mut q = Self::one();
        for _ in 0..count {
            q = q * &a;
            a = a.square();
        }
        q
    }

    /// Compute the square root 'x' of the field element such that x*x = self
    ///
    /// This function is not constant time, as it uses the Tonelli-Shanks algorithm
    pub fn sqrt(&self) -> CtOption<Self> {
        // use Tonelli-Shanks algorithm with hardcoded values related to the prime p,
        // Z=11, S=96, Q=340282366920938463463374607431768211455 (0xffffffffffffffffffffffffffffffff)
        // as such that p-1 is q*2^s and z the first non-quadratic residue in p.
        //
        // also (Q+1)/2 = 0x80000000000000000000000000000000 (1 bit to 1, 127 bits to 0)

        if !bool::from(self.is_quadratic_residue()) {
            return CtOption::from((CtEqual::ct_ne(self, self), FieldElement::zero()));
        }

        let n = self;
        let s = 96;

        let z = FieldElement::from_u64(11);

        let mut m = s;
        let mut c = z.square_add_rep(128); // c = z.power(&q)
        let mut t = n.square_add_rep(128); // t = n.power(&q)
        let mut r = n.square_rep(127); // r = n.power((q+1)/2)

        let one = FieldElement::one();

        while t != one {
            let mut tt = t.clone();
            let mut i = 0;
            while tt != one {
                tt = tt.square();
                i += 1;
                if i == m {
                    return CtOption::from((CtEqual::ct_ne(self, self), one));
                }
            }
            let e = m - i - 1;
            let b = if e == 0 { c.clone() } else { c.square_rep(e) };

            let b2 = b.square();
            r = r * b;
            t = t * &b2;
            c = b2;
            m = i;
        }

        let r2 = &r * &r;
        CtOption::from((CtEqual::ct_eq(&r2, self), r))
    }
}

fiat_field_ops_impl!(
    #[doc = "Element of the prime field Fp for scalar where p is the order of the SECP224R1 curve"]
    Scalar,
    224,
    ORDER_LIMBS,
    GM_LIMBS_SIZE,
    fiat_p224r1_scalar_nonzero,
    fiat_p224r1_scalar_add,
    fiat_p224r1_scalar_sub,
    fiat_p224r1_scalar_mul,
    fiat_p224r1_scalar_square,
    fiat_p224r1_scalar_opp,
    fiat_p224r1_scalar_to_bytes,
    fiat_p224r1_scalar_from_bytes,
    montgomery {
        fiat_p224r1_scalar_to_montgomery,
        fiat_p224r1_scalar_from_montgomery
    }
);

impl Scalar {
    /// Get the multiplicative inverse
    ///
    /// Note that 0 doesn't have a multiplicative inverse and will result in a panic
    /// TODO this will change to being a method of NonZeroScalar
    pub fn inverse(&self) -> Self {
        assert!(!self.is_zero());
        // 1*112,0*3,1*1,0*1,1*2,0*1,1*1,0*1,1*1,0*3,1*1,0*1,1*3,0*5,1*1,0*1,1*3,0*3,1*4,0*6,1*5,0*4,1*1,0*2,1*4,0*1,1*3,0*1,1*1,0*2,1*1,0*1,1*1,0*2,1*1,0*1,1*1,0*3,1*1,0*1,1*1,0*1,1*1,0*1,1*3,0*3,1*1,0*1,1*3,0*4,1*1,0*1,1*1,0*1,1*1,0*3,1*3,0*1,1*2
        let x2 = self.square() * self;
        let x3 = x2.square() * self;
        let x4 = x3.square() * self;
        let x5 = x4.square() * self;
        let x6 = x5.square() * self;
        let x9 = x6.square_rep(3) * &x3;
        let x18 = x9.square_rep(9) * &x9;
        let x27 = x18.square_rep(9) * &x9;
        let x54 = x27.square_rep(27) * &x27;
        let x108 = x54.square_rep(54) * &x54;
        let x112 = x108.square_rep(4) * &x4;

        let mut t1 = x112.square_rep(4) * self; // 1*112 0*3 1*1
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(4) * self; // 0*3 1*1
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(6) * self; // 0*5 1*1
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(7) * &x4; // 0*3 1*4
        t1 = t1.square_rep(11) * &x5; // 0*6 1*5
        t1 = t1.square_rep(5) * self; // 0*4 1*1
        t1 = t1.square_rep(6) * &x4; // 0*2 1*4
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(3) * self; // 0*2 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(3) * self; // 0*2 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(4) * self; // 0*3 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(4) * self; // 0*3 1*1
        t1 = t1.square_rep(4) * &x3; // 0*1 1*3
        t1 = t1.square_rep(5) * self; // 0*4 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(2) * self; // 0*1 1*1
        t1 = t1.square_rep(6) * &x3; // 0*3 1*3
        t1 = t1.square_rep(3) * &x2; // 0*1 1*2
        t1
    }
}

fiat_define_weierstrass_curve!(FieldElement);
fiat_define_weierstrass_points!(FieldElement);

impl Point {
    fn add_or_double<'b>(&self, other: &'b Point) -> Point {
        Point(self.0.add_or_double(&other.0, Curve))
    }
    fn scale<'b>(&self, other: &'b Scalar) -> Self {
        Point(self.0.scale(&other.to_bytes(), Curve))
    }
}

#[cfg(test)]
mod tests {
    mod fe {
        use super::super::FieldElement;
        use crate::{fiat_field_sqrt_unittest, fiat_field_unittest};

        fiat_field_unittest!(FieldElement);
        fiat_field_sqrt_unittest!(FieldElement);
    }
    mod gm {
        use super::super::Scalar;
        use crate::fiat_field_unittest;
        fiat_field_unittest!(Scalar);
    }
}
