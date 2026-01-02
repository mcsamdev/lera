use super::core::Complex;
use crate::scalar::Float;
use std::iter::{Product, Sum};

impl<T: Float> Sum for Complex<T> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::ZERO, |acc, x| acc + x)
    }
}

impl<T: Float> Product for Complex<T> {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::ONE, |acc, x| acc * x)
    }
}

impl<T: Float> Default for Complex<T> {
    fn default() -> Self {
        Self::ZERO
    }
}

impl From<f32> for Complex<f32> {
    fn from(f: f32) -> Self {
        Self::new(f, 0.0)
    }
}

impl From<f64> for Complex<f64> {
    fn from(f: f64) -> Self {
        Self::new(f, 0.0)
    }
}

impl From<i32> for Complex<f32> {
    fn from(i: i32) -> Self {
        Self::new(i as f32, 0.0)
    }
}

impl From<i64> for Complex<f64> {
    fn from(i: i64) -> Self {
        Self::new(i as f64, 0.0)
    }
}

impl From<u32> for Complex<f32> {
    fn from(i: u32) -> Self {
        Self::new(i as f32, 0.0)
    }
}

impl From<u64> for Complex<f64> {
    fn from(i: u64) -> Self {
        Self::new(i as f64, 0.0)
    }
}
