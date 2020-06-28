use num_bigint::{BigInt, BigUint, ToBigInt};
use num_traits::identities::{One, Zero};
use std::convert::TryInto;

pub fn mod_inverse(a: &BigUint, n: &BigUint) -> BigUint {
    let mut t = BigInt::zero();
    let mut newt = BigInt::one();
    let mut r: BigInt = n.to_bigint().unwrap();
    let mut newr: BigInt = a.to_bigint().unwrap();

    while !newr.is_zero() {
        let quotient = &r / &newr;

        let copyt = newt.clone();
        newt = &t - &quotient * newt;
        t = copyt;

        let copyr = newr.clone();
        newr = &r - &quotient * newr;
        r = copyr;
    }

    if r > BigInt::one() {
        panic!("a is not invertible")
    } else if t < BigInt::zero() {
        t = t + n.to_bigint().unwrap();
    }
    t.try_into().expect("no failure")
}
