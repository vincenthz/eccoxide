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

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum LegendreSymbol {
    One,
    Zero,
    MinusOne,
}

pub fn legendre_symbol(a: &BigUint, p: &BigUint) -> LegendreSymbol {
    use num_traits::cast::ToPrimitive;

    let r: u64 = a
        .modpow(&((p - BigUint::one()) / BigUint::from(2u64)), p)
        .to_u64()
        .unwrap();
    if r == 1 {
        LegendreSymbol::One
    } else if r == 0 {
        LegendreSymbol::Zero
    } else {
        LegendreSymbol::MinusOne
    }
}

// p need to be prime, but this is not checked
pub fn tonelli_shanks(n: &BigUint, p: &BigUint) -> Option<BigUint> {
    if legendre_symbol(n, p) != LegendreSymbol::One {
        return None;
    }

    let pm1 = p - BigUint::one();

    let mut q = pm1.clone();
    let mut s = 0;
    while (&q & BigUint::one()).is_zero() {
        s += 1;
        q >>= 1
    }

    // p == 3 mod 4 => obvious solution of n^((p+1)/4)
    if s == 1 {
        let r1 = n.modpow(&((p + BigUint::one()) / BigUint::from(4u32)), p);
        return Some(r1);
    }

    let mut z = BigUint::from(2u32);
    while z.modpow(&(&pm1 / BigUint::from(2u32)), p) != pm1 {
        z += BigUint::one()
    }
    let mut c = z.modpow(&q, p);

    let mut r = n.modpow(&((&q + BigUint::one()) / BigUint::from(2u32)), p);
    let mut t = n.modpow(&q, p);
    let mut m = s;

    while !t.is_one() {
        let mut tt = t.clone();
        let mut i = 0;
        while !tt.is_one() {
            tt = (&tt * &tt) % p;
            i += 1;
            if i == m {
                return None;
            }
        }
        let exponent = BigUint::from(2u32).modpow(&BigUint::from(m - i - 1u64), &pm1);
        let b = c.modpow(&exponent, p);
        let b2 = (&b * &b) % p;
        r = (r * b) % p;
        t = (t * &b2) % p;
        c = b2;
        m = i;
    }

    Some(r)
}
