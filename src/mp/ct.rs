//! Constant time operations
//!
//! This module exports traits to do basic checking operation in constant time,
//! those operations are:
//!
//! * CtZero : constant time zero and non-zero checking
//! * CtEqual : constant time equality and non-equality checking
//! * CtLesser : constant time less (<) and opposite greater-equal (>=) checking
//! * CtGreater : constant time greater (>) and opposite lesser-equal (<=) checking
//! * CtSelect : constant time selection between / conditional assignment of values
//!
//! And simple types to manipulate those capabilities in a safer way:
//!
//! * Choice : Constant time boolean and safe methods.
//!            this was initially called CtBool but aligned to other implementation.
//! * CtOption : Constant time Option type.
//!
//! Great care has been done to make operation constant so that it's useful in
//! cryptographic context, but we're not protected from implementation bug,
//! compiler optimisations, gamma rays and other Moon-Mars alignments.
//!
//! The general functionality would be a great addition to the rust core library
//! to have those type of things built-in and crucially more eyeballs.

/// Constant time boolean
///
/// This implementation uses a u64 under the hood, but it's never exposed
/// and only used through abstraction that push toward more constant time
/// operations.
///
/// Choice can be combined with simple And operation.
///
/// Choice can be converted back to a boolean operations, although
/// once this is done, the operation will likely be non-constant.
#[derive(Clone, Copy)]
pub struct Choice(pub(crate) u64);

/// Constant time equivalent to Option.
///
/// The T type is always present in the data structure,
/// it's just marked as valid / invalid with a Choice
/// type.
#[derive(Clone)]
pub struct CtOption<T> {
    present: Choice, // if present the value is there and valid
    t: T,
}

impl Choice {
    pub const TRUE: Self = Self(1);
    pub const FALSE: Self = Self(0);

    pub const fn is_true(self) -> bool {
        self.0 == 1
    }
    pub const fn is_false(self) -> bool {
        self.0 == 0
    }
    pub const fn negate(self) -> Self {
        Choice(1 ^ self.0)
    }

    // Return the choice as a `0`/`1` bit, in the shape expected by fiat-crypto
    // `*_selectznz` / `*_cmovznz` routines which use a u1
    pub(crate) const fn to_u1(self) -> u8 {
        self.0 as u8
    }

    // Return an all-ones mask (`!0`) when true, or an all-zeros mask (`0`) when
    // false. Used to implement branch-free selection.
    pub(crate) const fn to_mask(self) -> u64 {
        0u64.wrapping_sub(self.0)
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

impl core::ops::BitOr for Choice {
    type Output = Choice;
    fn bitor(self, b: Choice) -> Choice {
        Choice(self.0 | b.0)
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

/// Check in constant time if the object is zero or non-zero
///
/// Note that zero means 0 with integer primitive, or for array of integer
/// it means all elements are 0
pub trait CtZero {
    fn ct_zero(&self) -> Choice;
    fn ct_nonzero(&self) -> Choice;
}

/// Check in constant time if the left object is greater than right object
///
/// This equivalent to the > operator found in the core library.
#[allow(unused)]
pub trait CtGreater: Sized {
    fn ct_gt(a: Self, b: Self) -> Choice;
    fn ct_le(a: Self, b: Self) -> Choice {
        Self::ct_gt(b, a)
    }
}

/// Check in constant time if the left object is lesser than right object
///
/// This equivalent to the < operator found in the core library.
pub trait CtLesser: Sized {
    fn ct_lt(a: Self, b: Self) -> Choice;
    #[allow(unused)]
    fn ct_ge(a: Self, b: Self) -> Choice {
        Self::ct_lt(b, a)
    }
}

/// Check in constant time if the left object is equal to the right object
///
/// This equivalent to the == operator found in the core library.
pub trait CtEqual<Rhs: ?Sized = Self> {
    fn ct_eq(&self, b: &Rhs) -> Choice;
    fn ct_ne(&self, b: &Rhs) -> Choice {
        self.ct_eq(b).negate()
    }
}

/// Select between two values, or conditionally assign one, in constant time.
///
/// [`ct_select`](Self::ct_select)`(cond, a, b)` returns `a` when `cond` is true
/// and `b` otherwise; [`ct_assign`](Self::ct_assign)`(cond, other)` overwrites
/// `self` with `other` when `cond` is true and leaves it unchanged otherwise.
/// Neither must branch on `cond`, so both are safe to drive with a
/// secret-dependent choice (e.g. a scalar bit).
pub trait CtSelect: Sized {
    /// Conditionally assign `other` to `self`: `self` becomes `other` when
    /// `cond` is true, and is left unchanged otherwise.
    fn ct_assign(&mut self, cond: Choice, other: &Self);

    /// Returns `a` when `cond` is true and `b` otherwise. This is the equivalent
    /// of the `if cond { a } else { b }` expression of the core library.
    fn ct_select(cond: Choice, a: &Self, b: &Self) -> Self;
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

impl<const N: usize> CtEqual for [u8; N] {
    fn ct_eq(&self, b: &[u8; N]) -> Choice {
        let mut acc = 0u64;
        for (x, y) in self.iter().zip(b.iter()) {
            acc |= (*x as u64) ^ (*y as u64);
        }
        acc.ct_zero()
    }
}
impl<const N: usize> CtEqual for [u64; N] {
    fn ct_eq(&self, b: &[u64; N]) -> Choice {
        let mut acc = 0u64;
        for (x, y) in self.iter().zip(b.iter()) {
            acc |= x ^ y;
        }
        acc.ct_zero()
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

impl CtSelect for u64 {
    fn ct_assign(&mut self, cond: Choice, other: &Self) {
        // when cond is true the mask is !0, so `self` is XORed with the bits
        // that differ from `other` (becoming `other`); when false the mask is 0
        // and `self` is left unchanged.
        let mask = cond.to_mask();
        *self ^= mask & (*self ^ *other);
    }

    fn ct_select(cond: Choice, a: &Self, b: &Self) -> Self {
        // mask is !0 when cond is true and 0 when false, so this is `b` with the
        // differing bits of `a` flipped in exactly when selecting `a`.
        let mask = cond.to_mask();
        b ^ (mask & (a ^ b))
    }
}

impl<const N: usize> CtSelect for [u64; N] {
    fn ct_select(cond: Choice, a: &Self, b: &Self) -> Self {
        let mut out = [0u64; N];
        for (o, (x, y)) in out.iter_mut().zip(a.iter().zip(b.iter())) {
            *o = u64::ct_select(cond, x, y);
        }
        out
    }

    fn ct_assign(&mut self, cond: Choice, other: &Self) {
        for (s, o) in self.iter_mut().zip(other.iter()) {
            s.ct_assign(cond, o);
        }
    }
}

impl<const N: usize> CtSelect for [u8; N] {
    fn ct_select(cond: Choice, a: &Self, b: &Self) -> Self {
        let mask = cond.to_mask() as u8;
        let mut out = [0u8; N];
        for (o, (x, y)) in out.iter_mut().zip(a.iter().zip(b.iter())) {
            *o = y ^ (mask & (x ^ y));
        }
        out
    }

    fn ct_assign(&mut self, cond: Choice, other: &Self) {
        let mask = cond.to_mask() as u8;
        for (s, o) in self.iter_mut().zip(other.iter()) {
            *s ^= mask & (*s ^ *o);
        }
    }
}

// big endian representation of a number, but also leading byte of a array being the MSB.
impl<const N: usize> CtLesser for &[u8; N] {
    fn ct_lt(a: Self, b: Self) -> Choice {
        let mut borrow = 0u8;
        for (x, y) in a.iter().rev().zip(b.iter().rev()) {
            let x1: i16 = ((*x as i16) - (borrow as i16)) - (*y as i16);
            let x2: i8 = (x1 >> 8) as i8;
            borrow = (0x0 - x2) as u8;
        }
        let borrow = borrow as u64;
        Choice((borrow | borrow.wrapping_neg()) >> 63)
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

    #[test]
    fn test_ct_less() {
        let a: [u8; 4] = [0u8, 1, 2, 3];
        assert_eq!(<&[u8; 4]>::ct_lt(&a, &[1, 1, 2, 3]).is_true(), true);
    }

    #[test]
    fn ct_select_u64() {
        assert_eq!(u64::ct_select(Choice::TRUE, &0xdead, &0xbeef), 0xdead);
        assert_eq!(u64::ct_select(Choice::FALSE, &0xdead, &0xbeef), 0xbeef);
        // extreme bit patterns
        assert_eq!(u64::ct_select(Choice::TRUE, &!0, &0), !0);
        assert_eq!(u64::ct_select(Choice::FALSE, &!0, &0), 0);
    }

    #[test]
    fn ct_select_arrays() {
        let a8: [u8; 4] = [1, 2, 3, 4];
        let b8: [u8; 4] = [5, 6, 7, 8];
        assert_eq!(<[u8; 4]>::ct_select(Choice::TRUE, &a8, &b8), a8);
        assert_eq!(<[u8; 4]>::ct_select(Choice::FALSE, &a8, &b8), b8);

        let a64: [u64; 3] = [10, 20, 30];
        let b64: [u64; 3] = [40, 50, 60];
        assert_eq!(<[u64; 3]>::ct_select(Choice::TRUE, &a64, &b64), a64);
        assert_eq!(<[u64; 3]>::ct_select(Choice::FALSE, &a64, &b64), b64);
    }

    #[test]
    fn ct_assign_u64() {
        let mut x = 0xbeefu64;
        x.ct_assign(Choice::FALSE, &0xdead);
        assert_eq!(x, 0xbeef, "false must leave self unchanged");
        x.ct_assign(Choice::TRUE, &0xdead);
        assert_eq!(x, 0xdead, "true must overwrite self");

        // array variant via the default method
        let mut a: [u64; 3] = [1, 2, 3];
        a.ct_assign(Choice::TRUE, &[4, 5, 6]);
        assert_eq!(a, [4, 5, 6]);
        a.ct_assign(Choice::FALSE, &[7, 8, 9]);
        assert_eq!(a, [4, 5, 6]);
    }
}
