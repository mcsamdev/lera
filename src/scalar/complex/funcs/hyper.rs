use crate::scalar::complex::Complex;
use crate::scalar::Float;

impl<T: Float> Complex<T> {
    /// Computes the complex hyperbolic sine.
    ///
    /// # Algorithm
    /// Uses the identity:
    /// `sinh(a + bi) = sinh(a)cos(b) + i*cosh(a)sin(b)`
    #[inline]
    pub fn sinh(self) -> Self {
        Self::new(
            self.real.sinh() * self.imaginary.cos(),
            self.real.cosh() * self.imaginary.sin(),
        )
    }

    /// Computes the complex hyperbolic cosine.
    ///
    /// # Algorithm
    /// Uses the identity:
    /// `cosh(a + bi) = cosh(a)cos(b) + i*sinh(a)sin(b)`
    #[inline]
    pub fn cosh(self) -> Self {
        Self::new(
            self.real.cosh() * self.imaginary.cos(),
            self.real.sinh() * self.imaginary.sin(),
        )
    }

    /// Computes the complex hyperbolic tangent.
    ///
    /// # Algorithm
    /// Uses the formula: `tanh(z) = sinh(z) / cosh(z)`
    /// Computed directly to avoid separate sinh/cosh calls.
    #[inline]
    pub fn tanh(self) -> Self {
        // tanh(z) = sinh(z) / cosh(z)
        // More numerically stable: tanh(a+bi) = (sinh(2a) + i*sin(2b)) / (cosh(2a) + cos(2b))
        let two_a = T::TWO * self.real;
        let two_b = T::TWO * self.imaginary;
        let denom = two_a.cosh() + two_b.cos();
        Self::new(two_a.sinh() / denom, two_b.sin() / denom)
    }

    /// Computes the complex inverse hyperbolic sine.
    ///
    /// # Algorithm
    /// Uses the identity:
    /// `asinh(z) = ln(z + √(z² + 1))`
    #[inline]
    pub fn asinh(self) -> Self {
        (self + (self * self + Self::ONE).sqrt()).ln()
    }

    /// Computes the complex inverse hyperbolic cosine.
    ///
    /// # Algorithm
    /// Uses the identity:
    /// `acosh(z) = ln(z + √(z - 1) * √(z + 1))`
    #[inline]
    pub fn acosh(self) -> Self {
        (self + (self - Self::ONE).sqrt() * (self + Self::ONE).sqrt()).ln()
    }

    /// Computes the complex inverse hyperbolic tangent.
    ///
    /// # Algorithm
    /// Uses the identity:
    /// `atanh(z) = (1/2) * ln((1 + z) / (1 - z))`
    ///
    /// # Singularities
    /// Undefined at `z = ±1`.
    #[inline]
    pub fn atanh(self) -> Self {
        ((Self::ONE + self).ln() - (Self::ONE - self).ln()) * T::ONE_HALF
    }
}
