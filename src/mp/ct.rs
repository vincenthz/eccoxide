pub type Borrow = u8;
pub type IBorrow = i8;

#[derive(Clone, Copy)]
pub struct Choice(pub(crate) u64);

#[derive(Clone)]
pub struct CtOption<T> {
    present: Choice, // if present the value is there and valid
    t: T,
}

impl Choice {
    pub fn is_true(self) -> bool {
        self.0 == 1
    }
    pub fn is_false(self) -> bool {
        self.0 == 0
    }
    pub fn negate(self) -> Self {
        Choice(1 ^ self.0)
    }
}

impl From<Choice> for bool {
    fn from(c: Choice) -> bool {
        c.is_true()
    }
}

impl core::ops::BitAnd for Choice {
    type Output = Choice;
    fn bitand(self, b: Choice) -> Choice {
        Choice(self.0 & b.0)
    }
}

impl<T> From<(Choice, T)> for CtOption<T> {
    fn from(c: (Choice, T)) -> CtOption<T> {
        CtOption {
            present: c.0,
            t: c.1,
        }
    }
}

impl<T> CtOption<T> {
    pub fn into_option(self) -> Option<T> {
        if self.present.is_true() {
            Some(self.t)
        } else {
            None
        }
    }
}

pub trait CtZero {
    fn ct_zero(&self) -> Choice;
    fn ct_nonzero(&self) -> Choice;
}

pub trait CtGreater: Sized {
    fn ct_gt(a: Self, b: Self) -> Choice;
    fn ct_le(a: Self, b: Self) -> Choice {
        Self::ct_gt(b, a)
    }
}
pub trait CtLesser: Sized {
    fn ct_lt(a: Self, b: Self) -> Choice;
    fn ct_ge(a: Self, b: Self) -> Choice {
        Self::ct_lt(b, a)
    }
}

pub trait CtEqual<Rhs: ?Sized = Self> {
    fn ct_eq(&self, b: &Rhs) -> Choice;
    fn ct_ne(&self, b: &Rhs) -> Choice {
        self.ct_eq(b).negate()
    }
}

impl CtZero for u64 {
    fn ct_zero(&self) -> Choice {
        Choice(1 ^ ((self | self.wrapping_neg()) >> 63))
    }
    fn ct_nonzero(&self) -> Choice {
        Choice((self | self.wrapping_neg()) >> 63)
    }
}

impl CtEqual for u64 {
    fn ct_eq(&self, b: &Self) -> Choice {
        Self::ct_zero(&(self ^ b))
    }
    fn ct_ne(&self, b: &Self) -> Choice {
        Self::ct_nonzero(&(self ^ b))
    }
}

impl CtLesser for u64 {
    fn ct_lt(a: Self, b: Self) -> Choice {
        Choice((a ^ ((a ^ b) | ((a - b) ^ b))) >> 63)
    }
}

impl CtGreater for u64 {
    fn ct_gt(a: Self, b: Self) -> Choice {
        Self::ct_lt(b, a)
    }
}

impl<const N: usize> CtZero for [u8; N] {
    fn ct_zero(&self) -> Choice {
        let mut acc = 0u64;
        for b in self.iter() {
            acc |= *b as u64
        }
        acc.ct_zero()
    }
    fn ct_nonzero(&self) -> Choice {
        let mut acc = 0u64;
        for b in self.iter() {
            acc |= *b as u64
        }
        acc.ct_nonzero()
    }
}

impl<const N: usize> CtZero for [u64; N] {
    fn ct_zero(&self) -> Choice {
        let mut acc = 0u64;
        for b in self.iter() {
            acc |= b
        }
        acc.ct_zero()
    }
    fn ct_nonzero(&self) -> Choice {
        let mut acc = 0u64;
        for b in self.iter() {
            acc |= b
        }
        acc.ct_nonzero()
    }
}

impl CtZero for [u64] {
    fn ct_zero(&self) -> Choice {
        let mut acc = 0u64;
        for b in self.iter() {
            acc |= b
        }
        acc.ct_zero()
    }
    fn ct_nonzero(&self) -> Choice {
        let mut acc = 0u64;
        for b in self.iter() {
            acc |= b
        }
        acc.ct_nonzero()
    }
}

impl CtEqual for [u64] {
    fn ct_eq(&self, b: &[u64]) -> Choice {
        assert_eq!(self.len(), b.len());
        let mut acc = 0u64;
        for (x, y) in self.iter().zip(b.iter()) {
            acc |= x ^ y;
        }
        acc.ct_zero()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ct_zero() {
        assert_eq!(0u64.ct_zero().is_true(), true);
        assert_eq!(1u64.ct_zero().is_false(), true);
    }
}
