use crate::scalar::complex::Complex;
use crate::scalar::Float;

impl<T: Float> Complex<T> {
    /// Computes the complex sine.
    ///
    /// # Algorithm
    /// Uses the identity:
    /// `sin(a + bi) = sin(a)cosh(b) + i*cos(a)sinh(b)`
    ///
    /// # Examples
    /// ```
    /// use lera::scalar::Complex;
    ///
    /// let z = Complex::new(0.0_f64, 0.0);
    /// let result = z.sin();
    /// assert!(result.real().abs() < 1e-10);
    /// assert!(result.imaginary().abs() < 1e-10);
    /// ```
    #[inline]
    pub fn sin(self) -> Self {
        Self::new(
            self.real.sin() * self.imaginary.cosh(),
            self.real.cos() * self.imaginary.sinh(),
        )
    }

    /// Computes the complex cosine.
    ///
    /// # Algorithm
    /// Uses the identity:
    /// `cos(a + bi) = cos(a)cosh(b) - i*sin(a)sinh(b)`
    #[inline]
    pub fn cos(self) -> Self {
        Self::new(
            self.real.cos() * self.imaginary.cosh(),
            -self.real.sin() * self.imaginary.sinh(),
        )
    }

    /// Computes the complex tangent.
    ///
    /// # Algorithm
    /// Uses the identity: `tan(z) = sin(z) / cos(z)`
    /// Optimized to avoid computing sin and cos separately:
    /// - `denom = cos²(a) + sinh²(b)`
    /// - `re = sin(a)cos(a) / denom`
    /// - `im = sinh(b)cosh(b) / denom`
    ///
    /// # Singularities
    /// Has poles at `z = π/2 + nπ` for integer n.
    #[inline]
    pub fn tan(self) -> Self {
        let (sin_a, cos_a) = self.real.sin_cos();
        let sinh_b = self.imaginary.sinh();
        let cosh_b = self.imaginary.cosh();
        let denom = cos_a * cos_a + sinh_b * sinh_b;
        Self::new(sin_a * cos_a / denom, sinh_b * cosh_b / denom)
    }

    /// Computes the complex arcsine (inverse sine).
    ///
    /// # Algorithm
    /// Uses the identity:
    /// `asin(z) = -i * ln(iz + √(1 - z²))`
    ///
    /// # Branch Cuts
    /// Branch cuts along the real axis for |x| > 1.
    #[inline]
    pub fn asin(self) -> Self {
        let iz = self.mul_i();
        let root = (Self::one() - self * self).sqrt();
        (iz + root).ln().mul_i_neg()
    }

    /// Computes the complex arccosine (inverse cosine).
    ///
    /// # Algorithm
    /// Uses the identity:
    /// `acos(z) = -i * ln(z + √(z² - 1))`
    ///
    /// # Branch Cuts
    /// Branch cuts along the real axis for |x| > 1.
    #[inline]
    pub fn acos(self) -> Self {
        let root = (self * self - Self::one()).sqrt();
        (self + root).ln().mul_i_neg()
    }

    /// Computes the complex arctangent (inverse tangent).
    ///
    /// # Algorithm
    /// Uses the identity:
    /// `atan(z) = (i/2) * ln((1 - iz) / (1 + iz))`
    ///
    /// # Singularities
    /// Undefined at `z = ±i`.
    #[inline]
    pub fn atan(self) -> Self {
        let iz = self.mul_i();
        let one = Self::one();
        ((one - iz).ln() - (one + iz).ln()).mul_i() * T::ONE_HALF
    }
}
