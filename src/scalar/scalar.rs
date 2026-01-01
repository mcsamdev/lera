use super::{Complex, Float};

pub trait Scalar {
    /// The underlying real type used for bases and exponents
    type Base: Float;
    fn zero() -> Self;
    fn one() -> Self;
    fn abs(self) -> Self::Base;
    fn recip(self) -> Self;
    fn mul_add(self, a: Self, b: Self) -> Self;
    fn exp(self) -> Self;
    fn exp2(self) -> Self;
    fn exp_m1(self) -> Self;
    fn ln(self) -> Self;
    fn ln_1p(self) -> Self;
    fn log2(self) -> Self;
    fn log10(self) -> Self;
    fn log(self, base: Self::Base) -> Self;
    fn powi(self, n: i32) -> Self;
    fn powf(self, n: Self::Base) -> Self;
    fn sqrt(self) -> Self;
    fn cbrt(self) -> Self;
    fn sin(self) -> Self;
    fn cos(self) -> Self;
    fn tan(self) -> Self;
    fn asin(self) -> Self;
    fn acos(self) -> Self;
    fn atan(self) -> Self;
    fn sin_cos(self) -> (Self, Self)
    where
        Self: Sized;
    fn sinh(self) -> Self;
    fn cosh(self) -> Self;
    fn tanh(self) -> Self;
    fn asinh(self) -> Self;
    fn acosh(self) -> Self;
    fn atanh(self) -> Self;
    fn to_degrees(self) -> Self;
    fn to_radians(self) -> Self;
}

impl<T: Float> Scalar for T {
    type Base = T;
    fn zero() -> Self {
        Float::zero()
    }
    fn one() -> Self {
        Float::one()
    }
    fn abs(self) -> Self {
        Float::abs(self)
    }
    fn recip(self) -> Self {
        Float::recip(self)
    }
    fn mul_add(self, a: Self, b: Self) -> Self {
        Float::mul_add(self, a, b)
    }

    fn exp(self) -> Self {
        Float::exp(self)
    }
    fn exp2(self) -> Self {
        Float::exp2(self)
    }
    fn exp_m1(self) -> Self {
        Float::exp_m1(self)
    }
    fn ln(self) -> Self {
        Float::ln(self)
    }

    fn ln_1p(self) -> Self {
        Float::ln_1p(self)
    }
    fn log2(self) -> Self {
        Float::log2(self)
    }
    fn log10(self) -> Self {
        Float::log10(self)
    }
    fn log(self, base: Self::Base) -> Self {
        Self::log(self, base)
    }
    fn powi(self, n: i32) -> Self {
        Float::powi(self, n)
    }
    fn powf(self, n: Self::Base) -> Self {
        Float::powf(self, n)
    }
    fn sqrt(self) -> Self {
        Float::sqrt(self)
    }
    fn cbrt(self) -> Self {
        Float::cbrt(self)
    }
    fn sin(self) -> Self {
        Float::sin(self)
    }
    fn cos(self) -> Self {
        Float::cos(self)
    }
    fn tan(self) -> Self {
        Float::tan(self)
    }
    fn asin(self) -> Self {
        Float::asin(self)
    }
    fn acos(self) -> Self {
        Float::acos(self)
    }
    fn atan(self) -> Self {
        Float::atan(self)
    }
    fn sin_cos(self) -> (Self, Self) {
        Float::sin_cos(self)
    }
    fn sinh(self) -> Self {
        Float::sinh(self)
    }
    fn cosh(self) -> Self {
        Float::cosh(self)
    }
    fn tanh(self) -> Self {
        Float::tanh(self)
    }
    fn asinh(self) -> Self {
        Float::asinh(self)
    }
    fn acosh(self) -> Self {
        Float::acosh(self)
    }
    fn atanh(self) -> Self {
        Float::atanh(self)
    }
    fn to_degrees(self) -> Self {
        Float::to_degrees(self)
    }
    fn to_radians(self) -> Self {
        Float::to_radians(self)
    }
}

impl<T: Float> Scalar for Complex<T> {
    type Base = T;

    fn zero() -> Self {
        Self::zero()
    }

    fn one() -> Self {
        Self::one()
    }

    fn abs(self) -> T {
        Self::norm_sqr(self).sqrt()
    }

    fn recip(self) -> Self {
        Self::recip(self)
    }

    fn mul_add(self, a: Self, b: Self) -> Self {
        Self::mul_add(self, a, b)
    }

    fn exp(self) -> Self {
        Self::exp(self)
    }

    fn exp2(self) -> Self {
        Self::exp2(self)
    }

    fn exp_m1(self) -> Self {
        Self::exp_m1(self)
    }

    fn ln(self) -> Self {
        Self::ln(self)
    }

    fn ln_1p(self) -> Self {
        Self::ln_1p(self)
    }

    fn log2(self) -> Self {
        Self::log2(self)
    }

    fn log10(self) -> Self {
        Self::log10(self)
    }

    fn log(self, base: Self::Base) -> Self {
        Complex::log(self, base)
    }

    fn powi(self, n: i32) -> Self {
        Complex::powi(self, n)
    }
    fn powf(self, n: Self::Base) -> Self {
        Complex::powf(self, n)
    }
    fn sqrt(self) -> Self {
        Complex::sqrt(self)
    }
    fn cbrt(self) -> Self {
        Complex::cbrt(self)
    }
    fn sin(self) -> Self {
        Complex::sin(self)
    }
    fn cos(self) -> Self {
        Complex::cos(self)
    }
    fn tan(self) -> Self {
        Complex::tan(self)
    }
    fn asin(self) -> Self {
        Complex::asin(self)
    }
    fn acos(self) -> Self {
        Complex::acos(self)
    }
    fn atan(self) -> Self {
        Complex::atan(self)
    }
    fn sin_cos(self) -> (Self, Self) {
        (Complex::sin(self), Complex::cos(self))
    }
    fn sinh(self) -> Self {
        Complex::sinh(self)
    }
    fn cosh(self) -> Self {
        Complex::cosh(self)
    }
    fn tanh(self) -> Self {
        Complex::tanh(self)
    }
    fn asinh(self) -> Self {
        Complex::asinh(self)
    }
    fn acosh(self) -> Self {
        Complex::acosh(self)
    }
    fn atanh(self) -> Self {
        Complex::atanh(self)
    }
    fn to_degrees(self) -> Self {
        Complex::to_degrees(self)
    }
    fn to_radians(self) -> Self {
        Complex::to_radians(self)
    }
}
