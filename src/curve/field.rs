use crate::mp::ct::{CtEqual, CtOption};
use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sign {
    Positive,
    Negative,
}

// support trait for prime field
pub trait Field<Output = Self>:
    Sized
    + 'static
    + Send
    + Sync
    + Clone
    + PartialEq
    + Eq
    + CtEqual
    + fmt::Debug
    + fmt::Display
    + From<u64>
    + Add<Output = Output>
    + Sub<Output = Output>
    + Mul<Output = Output>
    + Neg<Output = Output>
    + for<'a> Add<&'a Self, Output = Output>
    + for<'a> Sub<&'a Self, Output = Output>
    + for<'a> Mul<&'a Self, Output = Output>
{
    fn zero() -> Output;
    fn is_zero(&self) -> bool;
    fn one() -> Output;
    fn double(&self) -> Output;

    fn inverse(&self) -> Output;
    fn sign(&self) -> Sign;

    fn square(&self) -> Output;
    fn cube(&self) -> Output;
}

pub trait FieldSqrt: Field {
    fn sqrt(&self) -> CtOption<Self>;
}
