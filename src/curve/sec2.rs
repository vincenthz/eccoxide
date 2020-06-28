macro_rules! prime_curve {
    ($m: ident) => {
        pub mod $m {
            use super::super::helper::mod_inverse;
            use crate::params::sec2::$m::*;
            use crate::scalar_impl;
            use lazy_static;
            use num_bigint::BigUint;

            lazy_static! {
                static ref P: BigUint = BigUint::from_bytes_be(&P_BYTES);
            }
            scalar_impl!(&*P);
        }
    };
}

prime_curve!(p112r1);
prime_curve!(p112r2);
prime_curve!(p128r1);
prime_curve!(p128r2);
prime_curve!(p160k1);
prime_curve!(p160r1);
prime_curve!(p160r2);
prime_curve!(p192k1);
prime_curve!(p192r1);
prime_curve!(p224k1);
prime_curve!(p224r1);
prime_curve!(p256k1);
prime_curve!(p256r1);
prime_curve!(p384r1);
prime_curve!(p521r1);
