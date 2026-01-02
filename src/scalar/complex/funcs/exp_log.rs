use crate::scalar::complex::Complex;
use crate::scalar::Float;

impl<T: Float> Complex<T> {
    /// Computes the complex exponential e^z.
    ///
    /// # Algorithm (Euler's Formula)
    /// For z = a + bi:
    /// `e^z = e^a * (cos(b) + i*sin(b))`
    ///
    /// # Examples
    /// ```
    /// use lera::scalar::{Complex, Float};
    ///
    /// // e^(iπ) = -1
    /// let z = Complex::new(0.0_f64, f64::PI);
    /// let result = z.exp();
    /// assert!((result.real() + 1.0).abs() < 1e-10);
    /// assert!(result.imaginary().abs() < 1e-10);
    /// ```
    #[inline]
    pub fn exp(self) -> Self {
        let ea = self.real.exp();
        let (s, c) = self.imaginary.sin_cos();
        Self::new(ea * c, ea * s)
    }

    /// Computes 2^z for complex z.
    ///
    /// # Algorithm
    /// Uses the identity: `2^z = e^(z * ln(2))`
    #[inline]
    pub fn exp2(self) -> Self {
        (self * T::LN_2).exp()
    }

    /// Computes e^z - 1.
    ///
    /// # Note
    /// For small |z|, this is more accurate than computing `exp(z) - 1`
    /// directly due to catastrophic cancellation.
    #[inline]
    pub fn exp_m1(self) -> Self {
        self.exp() - T::one()
    }

    /// Computes the principal value of the complex natural logarithm.
    ///
    /// # Algorithm
    /// For z = r * e^(iθ):
    /// `ln(z) = ln(r) + iθ`
    ///
    /// where r = |z| and θ = arg(z) ∈ (-π, π].
    ///
    /// # Branch Cut
    /// The branch cut is along the negative real axis.
    ///
    /// # Edge Cases
    /// - `ln(0)` returns `-∞ + 0i`
    /// - `ln(-1)` returns `0 + πi`
    ///
    /// # Examples
    /// ```
    /// use lera::scalar::{Complex, Float};
    ///
    /// let z = Complex::new(1.0_f64, 0.0);
    /// let result = z.ln();
    /// assert!(result.real().abs() < 1e-10); // ln(1) = 0
    /// ```
    #[inline]
    pub fn ln(self) -> Self {
        let r = self.real.hypot(self.imaginary);
        let theta = self.arg();
        Self::new(r.ln(), theta)
    }

    /// Computes ln(1 + z).
    ///
    /// # Note
    /// For small |z|, this may be more accurate than `(self + 1).ln()`
    /// due to precision loss in the addition.
    #[inline]
    pub fn ln_1p(self) -> Self {
        (self + T::one()).ln()
    }

    /// Computes the base-2 logarithm.
    ///
    /// # Algorithm
    /// Uses the change of base formula: `log₂(z) = ln(z) / ln(2)`
    #[inline]
    pub fn log2(self) -> Self {
        self.ln() / T::LN_2
    }

    /// Computes the base-10 logarithm.
    ///
    /// # Algorithm
    /// Uses the change of base formula: `log₁₀(z) = ln(z) / ln(10)`
    #[inline]
    pub fn log10(self) -> Self {
        self.ln() / T::LN_10
    }

    /// Computes the logarithm with a real base.
    ///
    /// # Algorithm
    /// Uses the change of base formula: `log_b(z) = ln(z) / ln(b)`
    ///
    /// # Arguments
    /// - `base`: The real logarithm base.
    #[inline]
    pub fn log(self, base: T) -> Self {
        self.ln() / base.ln()
    }

    /// Computes the logarithm with a complex base.
    ///
    /// # Algorithm
    /// Uses the change of base formula: `log_w(z) = ln(z) / ln(w)`
    ///
    /// # Arguments
    /// - `base`: The complex logarithm base.
    #[inline]
    pub fn logc(self, base: Self) -> Self {
        self.ln() / base.ln()
    }

    /// Raises the complex number to an integer power.
    ///
    /// # Algorithm
    /// Uses binary exponentiation (exponentiation by squaring) for O(log n) complexity.
    /// For negative exponents, computes the positive power then takes the reciprocal.
    ///
    /// # Examples
    /// ```
    /// use lera::scalar::Complex;
    ///
    /// let z = Complex::new(0.0_f64, 1.0); // i
    /// let result = z.powi(2);
    /// // i² = -1
    /// assert!((result.real() + 1.0).abs() < 1e-10);
    /// assert!(result.imaginary().abs() < 1e-10);
    /// ```
    #[inline]
    pub fn powi(self, exp: i32) -> Self {
        if exp == 0 {
            return Self::one();
        }
        let mut base = self;
        let mut e = exp.unsigned_abs();
        let mut result = Self::one();
        while e > 0 {
            if e & 1 == 1 {
                result *= base;
            }
            base = base * base;
            e >>= 1;
        }
        if exp < 0 { result.recip() } else { result }
    }

    /// Raises the complex number to a real power.
    ///
    /// # Algorithm
    /// Uses the identity: `z^n = e^(n * ln(z))`
    ///
    /// # Branch Cut
    /// Inherits the branch cut from [`ln`](Self::ln) (negative real axis).
    #[inline]
    pub fn powf(self, exp: T) -> Self {
        (self.ln() * exp).exp()
    }

    /// Raises the complex number to a complex power.
    ///
    /// # Algorithm
    /// Uses the identity: `z^w = e^(w * ln(z))`
    ///
    /// # Branch Cut
    /// Inherits the branch cut from [`ln`](Self::ln) (negative real axis).
    #[inline]
    pub fn powc(self, exp: Self) -> Self {
        (self.ln() * exp).exp()
    }
}
