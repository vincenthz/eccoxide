pub mod p256k1 {
    use super::super::helper::mod_inverse;
    use crate::params::sec2::p256k1::*;
    use crate::scalar_impl;
    use lazy_static;
    use num_bigint::BigUint;

    lazy_static! {
        static ref P: BigUint = BigUint::from_bytes_be(&P_BYTES);
    }
    scalar_impl!(&*P);
}
