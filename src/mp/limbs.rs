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
    pub fn iter_from_high(&self) -> core::iter::Rev<core::slice::Iter<'a, u64>> {
        self.0.iter().rev()
    }

    pub fn iter_from_low(&self) -> core::slice::Iter<'a, u64> {
        self.0.iter()
    }
}

impl<'a> LimbsBE<'a> {
    pub fn len(&self) -> usize {
        self.0.len()
    }
    pub fn iter_from_high(&self) -> core::slice::Iter<'a, u64> {
        self.0.iter()
    }

    pub fn iter_from_low(&self) -> core::iter::Rev<core::slice::Iter<'a, u64>> {
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

/// Pack big-endian bytes into saturated little-endian 64-bit limbs.
///
/// `L` must be wide enough for `N` bytes; any limb the bytes do not reach
/// stays zero.
pub(crate) const fn limbs_from_be<const N: usize, const L: usize>(bytes: &[u8; N]) -> [u64; L] {
    let mut out = [0u64; L];
    let mut i = 0;
    while i < N {
        // byte `i` counting up from the least significant end
        out[i / 8] |= (bytes[N - 1 - i] as u64) << ((i % 8) * 8);
        i += 1;
    }
    out
}

/// Pack little-endian bytes into saturated little-endian 64-bit limbs.
pub(crate) const fn limbs_from_le<const N: usize, const L: usize>(bytes: &[u8; N]) -> [u64; L] {
    let mut out = [0u64; L];
    let mut i = 0;
    while i < N {
        out[i / 8] |= (bytes[i] as u64) << ((i % 8) * 8);
        i += 1;
    }
    out
}

/// Unpack saturated little-endian 64-bit limbs into `N` little-endian bytes.
///
/// Bytes beyond `N` are dropped, so the caller must know the value fits — it
/// does when the value is reduced and `N` is the field's encoded size.
pub(crate) const fn limbs_to_le<const L: usize, const N: usize>(limbs: &[u64; L]) -> [u8; N] {
    let mut out = [0u8; N];
    let mut i = 0;
    while i < N {
        out[i] = (limbs[i / 8] >> ((i % 8) * 8)) as u8;
        i += 1;
    }
    out
}

/// `bytes - k` over a big-endian byte string
pub const fn be_sub_small<const N: usize>(bytes: &[u8; N], k: u8) -> [u8; N] {
    let mut out = *bytes;
    let mut borrow = k;
    let mut i = N;
    while i > 0 && borrow != 0 {
        i -= 1;
        if out[i] >= borrow {
            out[i] -= borrow;
            borrow = 0;
        } else {
            out[i] = (out[i] as u16 + 256 - borrow as u16) as u8;
            borrow = 1;
        }
    }
    out
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
    #[test]
    fn be_sub_small_borrows() {
        assert_eq!(be_sub_small(&[0x00, 0xab], 2), [0x00, 0xa9]);
        assert_eq!(be_sub_small(&[0x01, 0x00], 2), [0x00, 0xfe]);
        assert_eq!(
            be_sub_small(&[0x01, 0x00, 0x00, 0x01], 2),
            [0x00, 0xff, 0xff, 0xff]
        );
        assert_eq!(
            be_sub_small(&[0xff, 0x00, 0x00, 0x00], 1),
            [0xfe, 0xff, 0xff, 0xff]
        );
        assert_eq!(be_sub_small(&[0x00, 0x02], 2), [0x00, 0x00]);
    }
}
