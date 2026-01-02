use crate::scalar::complex::core::Complex;
use crate::scalar::Float;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

impl<T: Float> Neg for Complex<T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            real: -self.real,
            imaginary: -self.imaginary,
        }
    }
}

impl<T: Float> Add for Complex<T> {
    type Output = Self;
    fn add(self, other: Self) -> Self::Output {
        Self {
            real: self.real + other.real,
            imaginary: self.imaginary + other.imaginary,
        }
    }
}

impl<T: Float> Add<T> for Complex<T> {
    type Output = Self;
    fn add(self, other: T) -> Self::Output {
        Self {
            real: self.real + other,
            imaginary: self.imaginary,
        }
    }
}

impl<T: Float> AddAssign for Complex<T> {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl<T: Float> AddAssign<T> for Complex<T> {
    fn add_assign(&mut self, other: T) {
        *self = *self + other;
    }
}

impl<T: Float> Sub for Complex<T> {
    type Output = Self;
    fn sub(self, other: Self) -> Self::Output {
        Self {
            real: self.real - other.real,
            imaginary: self.imaginary - other.imaginary,
        }
    }
}

impl<T: Float> Sub<T> for Complex<T> {
    type Output = Self;
    fn sub(self, other: T) -> Self::Output {
        Self {
            real: self.real - other,
            imaginary: self.imaginary,
        }
    }
}

impl<T: Float> SubAssign for Complex<T> {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl<T: Float> SubAssign<T> for Complex<T> {
    fn sub_assign(&mut self, other: T) {
        *self = *self - other;
    }
}

impl<T: Float> Mul for Complex<T> {
    type Output = Self;
    fn mul(self, other: Self) -> Self::Output {
        Self {
            real: self.real * other.real - self.imaginary * other.imaginary,
            imaginary: self.real * other.imaginary + self.imaginary * other.real,
        }
    }
}

impl<T: Float> Mul<T> for Complex<T> {
    type Output = Self;
    fn mul(self, other: T) -> Self::Output {
        Self {
            real: self.real * other,
            imaginary: self.imaginary * other,
        }
    }
}

impl<T: Float> MulAssign for Complex<T> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<T: Float> MulAssign<T> for Complex<T> {
    fn mul_assign(&mut self, rhs: T) {
        self.real = self.real * rhs;
        self.imaginary = self.imaginary * rhs;
    }
}

impl<T: Float> Div for Complex<T> {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        let denom = rhs.norm_sqr();
        let num = self * rhs.conjugate();
        Self::new(num.real / denom, num.imaginary / denom)
    }
}

impl<T: Float> Div<T> for Complex<T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        Self::new(self.real / rhs, self.imaginary / rhs)
    }
}

impl<T: Float> DivAssign for Complex<T> {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl<T: Float> DivAssign<T> for Complex<T> {
    fn div_assign(&mut self, rhs: T) {
        self.real = self.real / rhs;
        self.imaginary = self.imaginary / rhs;
    }
}
