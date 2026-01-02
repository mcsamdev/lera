use crate::scalar::Complex;
use crate::scalar::Float;

impl<T: Float> Complex<T> {
    /// Computes the principal square root.
    ///
    /// # Algorithm
    /// Uses the formula for the principal square root:
    /// - `re = √((|z| + a) / 2)`
    /// - `im = b / (2 * re)` (with sign preserved)
    ///
    /// where z = a + bi.
    ///
    /// # Branch Cut
    /// The branch cut is along the negative real axis. The result always
    /// has non-negative real part.
    ///
    /// # Examples
    /// ```
    /// use lera::scalar::Complex;
    ///
    /// let z = Complex::new(-1.0_f64, 0.0);
    /// let result = z.sqrt();
    /// // √(-1) = i
    /// assert!(result.real().abs() < 1e-10);
    /// assert!((result.imaginary() - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn sqrt(self) -> Self {
        if self.real == T::zero() && self.imaginary == T::zero() {
            return self;
        }
        let r = self.real.hypot(self.imaginary);
        let two = T::one() + T::one();
        let re = ((r + self.real) / two).sqrt();
        let im = if re == T::zero() {
            r.sqrt()
        } else {
            self.imaginary / (two * re)
        };
        Self::new(re, im)
    }

    /// Computes the principal cube root.
    ///
    /// # Algorithm
    /// Converts to polar form (r, θ), computes `(r^(1/3), θ/3)`, then converts back.
    ///
    /// # Note
    /// Returns only the principal root. The other two cube roots can be obtained
    /// by multiplying by `e^(2πi/3)` and `e^(4πi/3)`.
    #[inline]
    pub fn cbrt(self) -> Self {
        if self.real == T::zero() && self.imaginary == T::zero() {
            return self;
        }
        let r = self.real.hypot(self.imaginary);
        let three = T::one() + T::one() + T::one();
        let theta_third = self.arg() / three;
        let r_cbrt = r.cbrt();
        Self::from_polar(r_cbrt, theta_third)
    }
}
