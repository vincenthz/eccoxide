macro_rules! prime_curve {
    ($m: ident, $szfe: expr) => {
        pub mod $m {
            use super::super::helper::mod_inverse;
            use crate::params::sec2::$m::*;
            use crate::{point_impl, scalar_impl};
            use lazy_static;
            use num_bigint::BigUint;

            lazy_static! {
                static ref P: BigUint = BigUint::from_bytes_be(&P_BYTES);
                static ref ORDER: BigUint = BigUint::from_bytes_be(&ORDER_BYTES);
            }
            scalar_impl!(&*P, $szfe);
            point_impl!(&*GX, &*GY);
        }
    };
}

prime_curve!(p112r1, 14);
prime_curve!(p112r2, 14);
prime_curve!(p128r1, 16);
prime_curve!(p128r2, 16);
prime_curve!(p160k1, 20);
prime_curve!(p160r1, 20);
prime_curve!(p160r2, 20);
prime_curve!(p192k1, 24);
prime_curve!(p192r1, 24);
prime_curve!(p224k1, 28);
prime_curve!(p224r1, 28);
prime_curve!(p256k1, 32);
prime_curve!(p256r1, 32);
prime_curve!(p384r1, 48);
prime_curve!(p521r1, 66);
