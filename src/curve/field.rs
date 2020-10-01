//! Field abstraction for abstracting basic field operations 0/1/+/*
//!
//! Also define some extended operations related like double (addition of self),
//! square (multiplication of self), and multiplicative inverse.
//!
//! We also expect a way to define a sign in the field, additive inverse (opposite), and the
//! overall subtraction operations so that the math checks out.
//!
//! For now the the field abstracted is expected to be a prime field (expected multiplicative inverse).
//! In the future, this will change to make the distinction
//! between field and prime field.
//!

use crate::mp::ct::{CtEqual, CtOption};
use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

/// Sign of a field element
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sign {
    Positive,
    Negative,
}

/// abstract trait for prime field support
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
