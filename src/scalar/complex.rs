// implementation of complex numbers for f32 and f64

/*
 * Type for complex numbers, the goal is to implement all possible or defined operations on complex numbers.
 * Examples of exclusions are:
 */
use super::Float;
use std::iter::{Product, Sum};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Complex<T: Float> {
    real: T,
    imaginary: T,
}

impl<T: Float> Complex<T> {
    pub fn new(real: T, imaginary: T) -> Self {
        Self { real, imaginary }
    }

    pub fn zero() -> Self {
        Self::new(T::zero(), T::zero())
    }

    pub fn one() -> Self {
        Self::new(T::one(), T::zero())
    }

    pub fn conjugate(self) -> Self {
        Self {
            real: self.real,
            imaginary: -self.imaginary,
        }
    }

    pub fn mul_i(self) -> Self {
        Self::new(-self.imaginary, self.real)
    }

    pub fn mul_i_neg(self) -> Self {
        Self::new(-self.imaginary, self.real)
    }

    pub fn magnitude(self) -> T {
        ((self.real * self.real) + (self.imaginary * self.imaginary)).sqrt()
    }

    pub fn recip(self) -> Self {
        let denom = Self::norm_sqr(self);
        Self::new(self.real / denom, -self.imaginary / denom)
    }

    pub fn mul_add(self, a: Self, b: Self) -> Self {
        let x = self.real;
        let y = self.imaginary;
        let u = a.real;
        let v = a.imaginary;

        // re_mul = x*u + (-y)*v
        let re_mul = x.mul_add(u, (-y) * v);

        // im_mul = x*v + y*u
        let im_mul = x.mul_add(v, y * u);

        Self::new(
            re_mul.mul_add(T::one(), b.real),
            im_mul.mul_add(T::one(), b.imaginary),
        )
    }

    pub fn mul_add_real(self, a: T, b: Self) -> Self {
        Self::new(
            self.real.mul_add(a, b.real),
            self.imaginary.mul_add(a, b.imaginary),
        )
    }

    pub fn arg(self) -> T {
        self.imaginary.atan2(self.real)
    }

    pub fn norm_sqr(self) -> T {
        self.real * self.real + self.imaginary * self.imaginary
    }

    pub fn from_polar(r: T, theta: T) -> Self {
        let (s, c) = theta.sin_cos();
        Self::new(r * c, r * s)
    }

    pub fn exp(self) -> Self {
        let ea = self.real.exp();
        let (s, c) = self.imaginary.sin_cos();
        Self::new(ea * c, ea * s)
    }

    pub fn exp2(self) -> Self {
        (self * T::LN_2).exp()
    }

    pub fn exp_m1(self) -> Self {
        self.exp() - T::one()
    }

    pub fn ln(self) -> Self {
        let r = self.real.hypot(self.imaginary);
        let theta = self.arg();
        Self::new(r.ln(), theta)
    }

    pub fn ln_1p(self) -> Self {
        self.ln() + T::one()
    }

    pub fn log2(self) -> Self {
        self.ln() / T::LN_2
    }

    pub fn log10(self) -> Self {
        self.ln() / T::LN_10
    }

    pub fn log(self, base: T) -> Self {
        self.ln() * (base.ln().recip())
    }

    pub fn logc(self, base: Self) -> Self {
        self.ln() * base.ln()
    }

    pub fn powi(self, exp: i32) -> Self {
        if exp == 0 {
            return Self::one();
        }
        let mut base = self;
        let mut e = exp.abs();
        let mut result = Self::one();
        while e > 0 {
            if e & 1 == 1 {
                result = result * base;
            }
            base = base * base;
            e >>= 1;
        }
        if exp < 0 { result.recip() } else { result }
    }

    pub fn powf(self, exp: T) -> Self {
        (self.ln() * exp).exp()
    }

    pub fn sqrt(self) -> Self {
        if self.real == T::zero() && self.imaginary == T::zero() {
            return self;
        }
        let r = self.real.hypot(self.imaginary);
        let re = ((r + self.real) * T::ONE_HALF).sqrt();
        //                               two without creating a constant
        let im = self.imaginary / (re * T::one() + T::one());
        Self::new(re, im)
    }

    pub fn cbrt(self) -> Self {
        if self.real == T::zero() && self.imaginary == T::zero() {
            return self;
        }
        let r = self.real.hypot(self.imaginary);
        let t = self.arg() / T::one() + T::one() + T::one();
        let (s, c) = t.sin_cos();
        Self::new(r.cbrt() * c, t.cbrt() * s)
    }

    pub fn sin(self) -> Self {
        Self::new(
            self.real.sin() * self.imaginary.cosh(),
            self.real.cos() * self.imaginary.sinh(),
        )
    }

    pub fn cos(self) -> Self {
        Self::new(
            self.real.cos() * self.imaginary.cosh(),
            -self.real.sin() * self.imaginary.sinh(),
        )
    }

    pub fn tan(self) -> Self {
        let denom =
            self.real.cos() * self.real.cos() + self.imaginary.sinh() * self.imaginary.sinh();
        Self::new(
            self.real.sin() * self.real.cos() / denom,
            self.imaginary.sinh() * self.imaginary.cosh() / denom,
        )
    }

    pub fn asin(self) -> Self {
        let iz = Self::new(-self.imaginary, self.real);
        let root = (Self::one() - self * self).sqrt();
        (iz + root).ln().mul_i_neg()
    }

    pub fn acos(self) -> Self {
        let root = (self * self - Self::one()).sqrt();
        (self + root).ln().mul_i_neg()
    }

    pub fn atan(self) -> Self {
        let iz = Self::new(-self.imaginary, self.real);
        let one = Self::one();
        ((one - iz).ln() - (one + iz).ln()).mul_i() * T::ONE_HALF
    }

    pub fn sinh(self) -> Self {
        Self::new(
            self.real.sinh() * self.imaginary.cos(),
            self.real.cosh() * self.imaginary.sin(),
        )
    }

    pub fn cosh(self) -> Self {
        Self::new(
            self.real.cosh() * self.imaginary.cos(),
            self.real.sinh() * self.imaginary.sin(),
        )
    }

    pub fn tanh(self) -> Self {
        let denom =
            self.real.cosh() * self.real.cosh() + self.imaginary.sin() * self.imaginary.sin();
        Self::new(
            self.real.sinh() * self.real.cosh() / denom,
            self.imaginary.sin() * self.imaginary.cos() / denom,
        )
    }

    pub fn asinh(self) -> Self {
        (self + (self * self + Self::one()).sqrt()).ln()
    }

    pub fn acosh(self) -> Self {
        let one = Self::one();
        (self + (self - one).sqrt() * (self + one).sqrt()).ln()
    }

    pub fn atanh(self) -> Self {
        let one = Self::one();
        ((one + self).ln() - (one - self).ln()) * (T::ONE_HALF)
    }

    pub fn to_degrees(self) -> Self {
        let r = self.magnitude();
        let theta = self.arg().to_degrees();
        Self::from_polar(r, theta)
    }

    pub fn to_radians(self) -> Self {
        let r = self.magnitude();
        let theta = self.arg().to_radians();
        Self::from_polar(r, theta)
    }

    pub fn real(self) -> T {
        self.real
    }
    pub fn imaginary(self) -> T {
        self.imaginary
    }
}

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

impl<T: Float> MulAssign<T> for Complex<T> {
    fn mul_assign(&mut self, rhs: T) {
        self.real = self.real * rhs;
        self.imaginary = self.imaginary * rhs;
    }
}

impl<T: Float> Div for Complex<T> {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        // (a + bi) / (c + di)
        // = (a + bi)(c - di) / (c² + d²)

        let denom = rhs.real * rhs.real + rhs.imaginary * rhs.imaginary;

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

impl<T: Float> Sum for Complex<T> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<T: Float> Product for Complex<T> {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
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
