use std::f32;
use std::f64;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, Sub, SubAssign};

pub trait Float:
    Copy
    + Clone
    + Neg<Output = Self>
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
    + Rem
    + PartialEq
    + PartialOrd
    + Default
{
    fn zero() -> Self;
    fn one() -> Self;
    fn abs(self) -> Self;
    fn recip(self) -> Self;
    fn mul_add(self, a: Self, b: Self) -> Self;
    fn exp(self) -> Self;
    fn exp2(self) -> Self;
    fn exp_m1(self) -> Self;
    fn ln(self) -> Self;
    fn ln_1p(self) -> Self;
    fn log2(self) -> Self;
    fn log10(self) -> Self;
    fn log(self, base: Self) -> Self;
    fn powi(self, n: i32) -> Self;
    fn powf(self, n: Self) -> Self;
    fn sqrt(self) -> Self;
    fn cbrt(self) -> Self;
    fn hypot(self, y: Self) -> Self;
    fn sin(self) -> Self;
    fn cos(self) -> Self;
    fn tan(self) -> Self;
    fn asin(self) -> Self;
    fn acos(self) -> Self;
    fn atan(self) -> Self;
    fn atan2(self, y: Self) -> Self;
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

    // Methods for special value checking
    fn is_nan(self) -> bool;
    fn is_infinite(self) -> bool;
    fn is_finite(self) -> bool;

    const EPSILON: Self;
    const MIN: Self;
    const MIN_POSITIVE: Self;
    const MAX: Self;
    const MANTISSA_DIGITS: usize;
    const RADIX: u32;
    const MIN_EXP: i32;
    const MAX_EXP: i32;
    const NEG_INFINITY: Self;
    const POSITIVE_INFINITY: Self;
    const NAN: Self;
    const PI: Self;
    const TAU: Self;
    const E: Self;
    const INFINITY: Self;
    const DEG_TO_RAD: Self;
    const RAD_TO_DEG: Self;

    // Fractional PI constants
    const FRAC_PI_2: Self; // π/2
    const FRAC_PI_3: Self; // π/3
    const FRAC_PI_4: Self; // π/4
    const FRAC_PI_6: Self; // π/6
    const FRAC_PI_8: Self; // π/8
    const FRAC_1_PI: Self; // 1/π
    const FRAC_2_PI: Self; // 2/π
    const FRAC_2_SQRT_PI: Self; // 2/√π

    // Logarithm constants
    const LN_2: Self; // ln(2)
    const LN_10: Self; // ln(10)
    const LOG2_E: Self; // log₂(e)
    const LOG2_10: Self; // log₂(10)
    const LOG10_E: Self; // log₁₀(e)
    const LOG10_2: Self; // log₁₀(2)

    // Root constants
    const SQRT_2: Self; // √2
    const SQRT_3: Self; // √3
    const FRAC_1_SQRT_2: Self; // 1/√2
    const FRAC_1_SQRT_3: Self; // 1/√3
    const CBRT_2: Self; // ∛2
    const CBRT_3: Self; // ∛3

    // Other useful constants
    const PHI: Self; // Golden ratio φ = (1 + √5)/2
    const ONE_HALF: Self; // 1/2
    const ONE: Self;
    const TWO: Self;
    const ZERO: Self;
}

impl Float for f32 {
    fn zero() -> Self {
        f32::ZERO
    }
    fn one() -> Self {
        f32::ONE
    }
    fn abs(self) -> Self {
        self.abs()
    }
    fn recip(self) -> Self {
        1.0 / self
    }
    fn mul_add(self, a: Self, b: Self) -> Self {
        self * a + b
    }
    fn exp(self) -> Self {
        self.exp()
    }
    fn exp2(self) -> Self {
        self.exp2()
    }
    fn exp_m1(self) -> Self {
        self.exp_m1()
    }
    fn ln(self) -> Self {
        self.ln()
    }
    fn ln_1p(self) -> Self {
        self.ln_1p()
    }
    fn log2(self) -> Self {
        self.log2()
    }
    fn log10(self) -> Self {
        self.log10()
    }
    fn log(self, base: Self) -> Self {
        self.log(base)
    }
    fn powi(self, n: i32) -> Self {
        self.powi(n)
    }
    fn powf(self, n: Self) -> Self {
        self.powf(n)
    }
    fn sqrt(self) -> Self {
        self.sqrt()
    }
    fn cbrt(self) -> Self {
        self.cbrt()
    }

    fn hypot(self, y: Self) -> Self {
        self.hypot(y)
    }

    fn sin(self) -> Self {
        self.sin()
    }
    fn cos(self) -> Self {
        self.cos()
    }
    fn tan(self) -> Self {
        self.tan()
    }
    fn asin(self) -> Self {
        self.asin()
    }
    fn acos(self) -> Self {
        self.acos()
    }
    fn atan(self) -> Self {
        self.atan()
    }

    fn atan2(self, y: Self) -> Self {
        self.atan2(y)
    }
    fn sin_cos(self) -> (Self, Self) {
        self.sin_cos()
    }
    fn sinh(self) -> Self {
        self.sinh()
    }
    fn cosh(self) -> Self {
        self.cosh()
    }
    fn tanh(self) -> Self {
        self.tanh()
    }
    fn asinh(self) -> Self {
        self.asinh()
    }
    fn acosh(self) -> Self {
        self.acosh()
    }
    fn atanh(self) -> Self {
        self.atanh()
    }
    fn to_degrees(self) -> Self {
        self.to_degrees()
    }
    fn to_radians(self) -> Self {
        self.to_radians()
    }

    fn is_nan(self) -> bool {
        self.is_nan()
    }
    fn is_infinite(self) -> bool {
        self.is_infinite()
    }
    fn is_finite(self) -> bool {
        self.is_finite()
    }

    const EPSILON: Self = f32::EPSILON;
    const MIN: Self = f32::MIN;
    const MIN_POSITIVE: Self = f32::MIN_POSITIVE;
    const MAX: Self = f32::MAX;
    const MANTISSA_DIGITS: usize = f32::MANTISSA_DIGITS as usize;
    const RADIX: u32 = f32::RADIX;
    const MIN_EXP: i32 = f32::MIN_EXP;
    const MAX_EXP: i32 = f32::MAX_EXP;
    const NEG_INFINITY: Self = f32::NEG_INFINITY;
    const POSITIVE_INFINITY: Self = f32::INFINITY;
    const NAN: Self = f32::NAN;
    const PI: Self = f32::consts::PI;
    const TAU: Self = f32::consts::TAU;
    const E: Self = f32::consts::E;
    const INFINITY: Self = f32::INFINITY;
    const DEG_TO_RAD: Self = f32::consts::PI / 180.0;
    const RAD_TO_DEG: Self = 180.0 / f32::consts::PI;

    // Fractional PI constants
    const FRAC_PI_2: Self = f32::consts::FRAC_PI_2;
    const FRAC_PI_3: Self = f32::consts::FRAC_PI_3;
    const FRAC_PI_4: Self = f32::consts::FRAC_PI_4;
    const FRAC_PI_6: Self = f32::consts::FRAC_PI_6;
    const FRAC_PI_8: Self = f32::consts::FRAC_PI_8;
    const FRAC_1_PI: Self = f32::consts::FRAC_1_PI;
    const FRAC_2_PI: Self = f32::consts::FRAC_2_PI;
    const FRAC_2_SQRT_PI: Self = f32::consts::FRAC_2_SQRT_PI;

    // Logarithm constants
    const LN_2: Self = f32::consts::LN_2;
    const LN_10: Self = f32::consts::LN_10;
    const LOG2_E: Self = f32::consts::LOG2_E;
    const LOG2_10: Self = f32::consts::LOG2_10;
    const LOG10_E: Self = f32::consts::LOG10_E;
    const LOG10_2: Self = f32::consts::LOG10_2;

    // Root constants
    const SQRT_2: Self = f32::consts::SQRT_2;
    const SQRT_3: Self = 1.732050807568877293527446341505872367_f32;
    const FRAC_1_SQRT_2: Self = f32::consts::FRAC_1_SQRT_2;
    const FRAC_1_SQRT_3: Self = 0.577350269189625764509148780501957456_f32;
    const CBRT_2: Self = 1.259921049894873164767210607278228350_f32;
    const CBRT_3: Self = 1.442249570307408382321638310780109588_f32;

    // Other useful constants
    const PHI: Self = 1.618033988749894848204586834365638118_f32;
    const ONE_HALF: Self = 0.5_f32;
    const ONE: Self = 1.0_f32;
    const TWO: Self = 2.0_f32;
    const ZERO: Self = 0.0_f32;
}
impl Float for f64 {
    fn zero() -> Self {
        0.0
    }
    fn one() -> Self {
        1.0
    }
    fn abs(self) -> Self {
        self.abs()
    }
    fn recip(self) -> Self {
        1.0 / self
    }
    fn mul_add(self, a: Self, b: Self) -> Self {
        self * a + b
    }
    fn exp(self) -> Self {
        self.exp()
    }
    fn exp2(self) -> Self {
        self.exp2()
    }
    fn exp_m1(self) -> Self {
        self.exp_m1()
    }
    fn ln(self) -> Self {
        self.ln()
    }
    fn ln_1p(self) -> Self {
        self.ln_1p()
    }
    fn log2(self) -> Self {
        self.log2()
    }
    fn log10(self) -> Self {
        self.log10()
    }
    fn log(self, base: Self) -> Self {
        self.log(base)
    }
    fn powi(self, n: i32) -> Self {
        self.powi(n)
    }
    fn powf(self, n: Self) -> Self {
        self.powf(n)
    }
    fn sqrt(self) -> Self {
        self.sqrt()
    }
    fn cbrt(self) -> Self {
        self.cbrt()
    }
    fn hypot(self, y: Self) -> Self {
        self.hypot(y)
    }
    fn sin(self) -> Self {
        self.sin()
    }
    fn cos(self) -> Self {
        self.cos()
    }
    fn tan(self) -> Self {
        self.tan()
    }
    fn asin(self) -> Self {
        self.asin()
    }
    fn acos(self) -> Self {
        self.acos()
    }
    fn atan(self) -> Self {
        self.atan()
    }

    fn atan2(self, y: Self) -> Self {
        self.atan2(y)
    }
    fn sin_cos(self) -> (Self, Self) {
        self.sin_cos()
    }
    fn sinh(self) -> Self {
        self.sinh()
    }
    fn cosh(self) -> Self {
        self.cosh()
    }
    fn tanh(self) -> Self {
        self.tanh()
    }
    fn asinh(self) -> Self {
        self.asinh()
    }
    fn acosh(self) -> Self {
        self.acosh()
    }
    fn atanh(self) -> Self {
        self.atanh()
    }
    fn to_degrees(self) -> Self {
        self.to_degrees()
    }
    fn to_radians(self) -> Self {
        self.to_radians()
    }

    fn is_nan(self) -> bool {
        self.is_nan()
    }
    fn is_infinite(self) -> bool {
        self.is_infinite()
    }
    fn is_finite(self) -> bool {
        self.is_finite()
    }

    const EPSILON: Self = f64::EPSILON;
    const MIN: Self = f64::MIN;
    const MIN_POSITIVE: Self = f64::MIN_POSITIVE;
    const MAX: Self = f64::MAX;
    const MANTISSA_DIGITS: usize = f64::MANTISSA_DIGITS as usize;
    const RADIX: u32 = f64::RADIX;
    const MIN_EXP: i32 = f64::MIN_EXP;
    const MAX_EXP: i32 = f64::MAX_EXP;
    const NEG_INFINITY: Self = f64::NEG_INFINITY;
    const POSITIVE_INFINITY: Self = f64::INFINITY;
    const NAN: Self = f64::NAN;
    const PI: Self = f64::consts::PI;
    const TAU: Self = f64::consts::TAU;
    const E: Self = f64::consts::E;
    const INFINITY: Self = f64::INFINITY;
    const DEG_TO_RAD: Self = f64::consts::PI / 180.0;
    const RAD_TO_DEG: Self = 180.0 / f64::consts::PI;
    // Fractional PI constants
    const FRAC_PI_2: Self = f64::consts::FRAC_PI_2;
    const FRAC_PI_3: Self = f64::consts::FRAC_PI_3;
    const FRAC_PI_4: Self = f64::consts::FRAC_PI_4;
    const FRAC_PI_6: Self = f64::consts::FRAC_PI_6;
    const FRAC_PI_8: Self = f64::consts::FRAC_PI_8;
    const FRAC_1_PI: Self = f64::consts::FRAC_1_PI;
    const FRAC_2_PI: Self = f64::consts::FRAC_2_PI;
    const FRAC_2_SQRT_PI: Self = f64::consts::FRAC_2_SQRT_PI;

    // Logarithm constants
    const LN_2: Self = f64::consts::LN_2;
    const LN_10: Self = f64::consts::LN_10;
    const LOG2_E: Self = f64::consts::LOG2_E;
    const LOG2_10: Self = f64::consts::LOG2_10;
    const LOG10_E: Self = f64::consts::LOG10_E;
    const LOG10_2: Self = f64::consts::LOG10_2;

    // Root constants
    const SQRT_2: Self = f64::consts::SQRT_2;
    const SQRT_3: Self = 1.732050807568877293527446341505872367_f64;
    const FRAC_1_SQRT_2: Self = f64::consts::FRAC_1_SQRT_2;
    const FRAC_1_SQRT_3: Self = 0.577350269189625764509148780501957456_f64;
    const CBRT_2: Self = 1.259921049894873164767210607278228350_f64;
    const CBRT_3: Self = 1.442249570307408382321638310780109588_f64;

    // Other useful constants
    const PHI: Self = 1.618033988749894848204586834365638118_f64;
    const ONE_HALF: Self = 0.5_f64;
    const ONE: Self = 1.0_f64;
    const TWO: Self = 2.0_f64;
    const ZERO: Self = 0.0_f64;
}
