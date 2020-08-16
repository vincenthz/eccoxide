pub type Borrow = u8;
pub type IBorrow = i8;

#[derive(Clone, Copy)]
pub struct Choice(pub(super) u64);

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

pub trait CtZero: Sized {
    fn ct_zero(a: Self) -> Choice;
    fn ct_nonzero(a: Self) -> Choice;
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

pub trait CtEqual: Sized {
    fn ct_eq(a: Self, b: Self) -> Choice;
    fn ct_ne(a: Self, b: Self) -> Choice {
        Self::ct_eq(a, b).negate()
    }
}

impl CtZero for u64 {
    fn ct_zero(a: Self) -> Choice {
        Choice(1 ^ ((a | a.wrapping_neg()) >> 63))
    }
    fn ct_nonzero(a: Self) -> Choice {
        Choice((a | a.wrapping_neg()) >> 63)
    }
}

impl CtEqual for u64 {
    fn ct_eq(a: Self, b: Self) -> Choice {
        Self::ct_zero(a ^ b)
    }
    fn ct_ne(a: Self, b: Self) -> Choice {
        Self::ct_nonzero(a ^ b)
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

impl CtZero for &[u64] {
    fn ct_zero(a: &[u64]) -> Choice {
        let mut acc = 0;
        for b in a.iter() {
            acc |= b
        }
        u64::ct_zero(acc)
    }
    fn ct_nonzero(a: Self) -> Choice {
        let mut acc = 0;
        for b in a.iter() {
            acc |= b
        }
        u64::ct_nonzero(acc)
    }
}

impl CtEqual for &[u64] {
    fn ct_eq(a: &[u64], b: &[u64]) -> Choice {
        assert_eq!(a.len(), b.len());
        let mut acc = 0;
        for (x, y) in a.iter().zip(b.iter()) {
            acc |= x ^ y;
        }
        u64::ct_zero(acc)
    }
}
