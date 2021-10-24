#![allow(dead_code)]

use super::ct::*;

pub type Limb = u64;

pub type Borrow = u8;
pub type IBorrow = i8;

pub struct LimbsLE<'a>(pub &'a [Limb]);

pub struct LimbsBE<'a>(pub &'a [Limb]);

impl<'a> LimbsLE<'a> {
    pub fn len(&self) -> usize {
        self.0.len()
    }
    pub fn iter_from_high(&self) -> core::iter::Rev<std::slice::Iter<'a, u64>> {
        self.0.iter().rev()
    }

    pub fn iter_from_low(&self) -> std::slice::Iter<'a, u64> {
        self.0.iter()
    }
}

impl<'a> LimbsBE<'a> {
    pub fn len(&self) -> usize {
        self.0.len()
    }
    pub fn iter_from_high(&self) -> std::slice::Iter<'a, u64> {
        self.0.iter()
    }

    pub fn iter_from_low(&self) -> core::iter::Rev<std::slice::Iter<'a, u64>> {
        self.0.iter().rev()
    }
}

// borrowed from fiat-crypto subborrow routine
fn limb_subborrow(out1: &mut u64, out2: &mut Borrow, arg1: Borrow, arg2: u64, arg3: u64) -> () {
    let x1: i128 = ((arg2 as i128) - (arg1 as i128)) - (arg3 as i128);
    let x2: IBorrow = (x1 >> 64) as IBorrow;
    let x3: u64 = (x1 & (0xffffffffffffffff as i128)) as u64;
    *out1 = x3;
    *out2 = ((0x0 as IBorrow) - (x2 as IBorrow)) as Borrow;
}

// Check that the value a is less than the value b
pub fn limbsbe_le<'a, 'b>(a: LimbsBE<'a>, b: LimbsBE<'b>) -> Choice {
    assert_eq!(a.len(), b.len());

    let mut borrow: Borrow = 0;
    let mut out = 0u64;
    for (x, y) in b.iter_from_low().zip(a.iter_from_low()) {
        let copied_borrow = borrow;
        limb_subborrow(&mut out, &mut borrow, copied_borrow, *x, *y);
    }
    Choice(1 ^ borrow as u64)
}

pub fn limbsbe_lt<'a, 'b>(a: LimbsBE<'a>, b: LimbsBE<'b>) -> Choice {
    assert_eq!(a.len(), b.len());

    let mut borrow: Borrow = 0;
    let mut out = 0u64;
    for (x, y) in a.iter_from_low().zip(b.iter_from_low()) {
        let copied_borrow = borrow;
        limb_subborrow(&mut out, &mut borrow, copied_borrow, *x, *y);
    }
    let borrow = borrow as u64;
    Choice((borrow | borrow.wrapping_neg()) >> 63)
}

pub fn limbsle_le<'a, 'b>(a: LimbsLE<'a>, b: LimbsLE<'b>) -> Choice {
    assert_eq!(a.len(), b.len());

    let mut borrow: Borrow = 0;
    let mut out = 0u64;
    for (x, y) in b.iter_from_low().zip(a.iter_from_low()) {
        let copied_borrow = borrow;
        limb_subborrow(&mut out, &mut borrow, copied_borrow, *x, *y)
    }
    Choice(1 ^ borrow as u64)
}

pub fn limbsle_lt<'a, 'b>(a: LimbsLE<'a>, b: LimbsLE<'b>) -> Choice {
    assert_eq!(a.len(), b.len());

    let mut borrow: Borrow = 0;
    let mut out = 0u64;
    for (x, y) in a.iter_from_low().zip(b.iter_from_low()) {
        let copied_borrow = borrow;
        limb_subborrow(&mut out, &mut borrow, copied_borrow, *x, *y);
    }
    let borrow = borrow as u64;
    Choice((borrow | borrow.wrapping_neg()) >> 63)
}

impl<'a> CtEqual for LimbsLE<'a> {
    fn ct_eq(&self, b: &Self) -> Choice {
        self.0.ct_eq(b.0)
    }
}

impl<'a> CtEqual for LimbsBE<'a> {
    fn ct_eq(&self, b: &Self) -> Choice {
        self.0.ct_eq(b.0)
    }
}

impl<'a> CtZero for LimbsLE<'a> {
    fn ct_zero(&self) -> Choice {
        self.0.ct_zero()
    }
    fn ct_nonzero(&self) -> Choice {
        self.0.ct_nonzero()
    }
}

impl<'a> CtLesser for LimbsLE<'a> {
    fn ct_lt(a: Self, b: Self) -> Choice {
        limbsle_lt(a, b)
    }
}

impl<'a> CtLesser for LimbsBE<'a> {
    fn ct_lt(a: Self, b: Self) -> Choice {
        limbsbe_lt(a, b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn le() {
        assert_eq!(
            false,
            limbsbe_le(
                LimbsBE(&[1, 0, 0, 0]),
                LimbsBE(&[0, u64::MAX, u64::MAX, u64::MAX])
            )
            .into()
        );
        assert_eq!(
            true,
            limbsbe_le(LimbsBE(&[1, 2, 3]), LimbsBE(&[1, 2, 3])).into()
        );
        assert_eq!(
            true,
            limbsbe_le(LimbsBE(&[1, 2, 3]), LimbsBE(&[1, 3, 3])).into()
        );
        assert_eq!(
            true,
            limbsbe_le(LimbsBE(&[0, 2, 3]), LimbsBE(&[1, 2, 3])).into(),
        );
        assert_eq!(
            false,
            limbsbe_le(LimbsBE(&[1, 4, 2]), LimbsBE(&[1, 2, 3])).into(),
        );
        assert_eq!(
            false,
            limbsbe_le(LimbsBE(&[2, 0, 2]), LimbsBE(&[1, 2, 3])).into(),
        );
    }

    #[test]
    fn lt() {
        assert_eq!(
            false,
            limbsbe_le(
                LimbsBE(&[1, 0, 0, 0]),
                LimbsBE(&[0, u64::MAX, u64::MAX, u64::MAX])
            )
            .into(),
        );
        assert_eq!(
            false,
            limbsbe_lt(LimbsBE(&[1, 2, 3]), LimbsBE(&[1, 2, 3])).into(),
        );
        assert_eq!(
            true,
            limbsbe_lt(LimbsBE(&[1, 2, 3]), LimbsBE(&[1, 3, 3])).into(),
        );
        assert_eq!(
            true,
            limbsbe_lt(LimbsBE(&[0, 2, 3]), LimbsBE(&[1, 2, 3])).into(),
        );
        assert_eq!(
            false,
            limbsbe_lt(LimbsBE(&[1, 4, 2]), LimbsBE(&[1, 2, 3])).into(),
        );
        assert_eq!(
            false,
            limbsbe_lt(LimbsBE(&[2, 0, 2]), LimbsBE(&[1, 2, 3])).into(),
        );
    }
}
