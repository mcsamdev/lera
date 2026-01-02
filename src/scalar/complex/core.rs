//! Implementation of complex numbers for f32 and f64.
//!
//! This module provides a generic [`crate::scalar::complex::complex::Complex<T> `] type that supports all standard
//! mathematical operations on complex numbers, where `T` is any type implementing
//! the [`Float`] trait.

use crate::scalar::Float;

/// A complex number represented as `real + imaginary * i`.
///
/// # Type Parameters
/// - `T`: The underlying floating-point type (must implement [`Float`]).
///
/// # Examples
/// ```
/// use lera::scalar::{Complex, Float};
///
/// let z = Complex::new(3.0_f64, 4.0);
/// assert_eq!(z.abs(), 5.0);
/// assert_eq!(z.conjugate(), Complex::new(3.0, -4.0));
/// ```
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Complex<T: Float> {
    pub(crate) real: T,
    pub(crate) imaginary: T,
}

impl<T: Float> Complex<T> {
    /// Creates a new complex number from real and imaginary parts.
    ///
    /// # Arguments
    /// - `real`: The real component.
    /// - `imaginary`: The imaginary component.
    ///
    /// # Examples
    /// ```
    /// use lera::scalar::Complex;
    ///
    /// let z = Complex::new(1.0_f64, 2.0);
    /// assert_eq!(z.real(), 1.0);
    /// assert_eq!(z.imaginary(), 2.0);
    /// ```
    #[inline]
    pub const fn new(real: T, imaginary: T) -> Self {
        Self { real, imaginary }
    }

    /// Returns the complex conjugate (a - bi for a + bi).
    ///
    /// The conjugate satisfies: `z * z.conjugate() == |z|²`
    ///
    /// # Examples
    /// ```
    /// use lera::scalar::Complex;
    ///
    /// let z = Complex::new(3.0_f64, 4.0);
    /// assert_eq!(z.conjugate(), Complex::new(3.0, -4.0));
    /// ```
    #[inline]
    pub fn conjugate(self) -> Self {
        Self {
            real: self.real,
            imaginary: -self.imaginary,
        }
    }

    /// Multiplies by the imaginary unit i.
    ///
    /// Equivalent to `self * i`, computed as `(-imaginary, real)`.
    ///
    /// # Mathematical Identity
    /// For z = a + bi: `z * i = -b + ai`
    #[inline]
    pub fn mul_i(self) -> Self {
        Self::new(-self.imaginary, self.real)
    }

    /// Multiplies by the negative imaginary unit -i.
    ///
    /// Equivalent to `self * (-i)`, computed as `(imaginary, -real)`.
    ///
    /// # Mathematical Identity
    /// For z = a + bi: `z * (-i) = b - ai`
    #[inline]
    pub fn mul_i_neg(self) -> Self {
        Self::new(self.imaginary, -self.real)
    }

    /// Returns the abs (absolute value) of the complex number.
    ///
    /// Computes `√(real² + imaginary²)` using the standard Euclidean norm.
    ///
    /// # Note
    /// For better numerical stability with extreme values, consider using
    /// [`norm_sqr`](Self::norm_sqr) when you only need the squared abs.
    ///
    /// # Examples
    /// ```
    /// use lera::scalar::Complex;
    ///
    /// let z = Complex::new(3.0_f64, 4.0);
    /// assert_eq!(z.abs(), 5.0);
    /// ```
    #[inline]
    pub fn abs(self) -> T {
        self.real.hypot(self.imaginary)
    }

    /// Returns the multiplicative inverse (1/z).
    ///
    /// # Algorithm
    /// Uses the formula: `1/(a+bi) = (a-bi)/(a²+b²)`
    ///
    /// # Panics
    /// May produce infinity or NaN for zero or near-zero inputs.
    ///
    /// # Examples
    /// ```
    /// use lera::scalar::Complex;
    ///
    /// let z = Complex::new(1.0_f64, 1.0);
    /// let inv = z.recip();
    /// // z * inv ≈ 1
    /// ```
    #[inline]
    pub fn recip(self) -> Self {
        let denom = self.norm_sqr();
        Self::new(self.real / denom, -self.imaginary / denom)
    }

    /// Computes `self * a + b` (fused multiply-add for complex numbers).
    ///
    /// # Algorithm
    /// Performs complex multiplication followed by addition:
    /// - Real part: `(x*u - y*v) + b.real`
    /// - Imaginary part: `(x*v + y*u) + b.imaginary`
    ///
    /// Uses fused multiply-add operations on the underlying floats for
    /// improved precision where available.
    #[inline]
    pub fn mul_add(self, a: Self, b: Self) -> Self {
        let x = self.real;
        let y = self.imaginary;
        let u = a.real;
        let v = a.imaginary;

        let re_mul = x.mul_add(u, (-y) * v);
        let im_mul = x.mul_add(v, y * u);

        Self::new(re_mul + b.real, im_mul + b.imaginary)
    }

    /// Computes `self * a + b` where `a` is a real scalar.
    ///
    /// More efficient than full complex multiply-add when the multiplier is real.
    #[inline]
    pub fn mul_add_real(self, a: T, b: Self) -> Self {
        Self::new(
            self.real.mul_add(a, b.real),
            self.imaginary.mul_add(a, b.imaginary),
        )
    }

    /// Returns the argument (phase angle) in radians.
    ///
    /// # Returns
    /// The angle θ in the range `(-π, π]` such that `z = |z| * e^(iθ)`.
    ///
    /// # Algorithm
    /// Uses `atan2(imaginary, real)` for correct quadrant handling.
    ///
    /// # Examples
    /// ```
    /// use lera::scalar::{Complex, Float};
    ///
    /// let z = Complex::new(1.0_f64, 1.0);
    /// assert!((z.arg() - f64::FRAC_PI_4).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn arg(self) -> T {
        self.imaginary.atan2(self.real)
    }

    /// Returns the squared abs (squared norm).
    ///
    /// Computes `real² + imaginary²` without taking the square root.
    ///
    /// # Use Cases
    /// - More efficient than `abs()` when you only need comparisons
    /// - Avoids potential overflow/underflow from the square root
    /// - Useful for computing `|z|²` directly
    #[inline]
    pub fn norm_sqr(self) -> T {
        let ar = self.real.abs();
        let ai = self.imaginary.abs();
        let m = if ar >= ai { ar } else { ai };
        if m == T::zero() {
            T::zero()
        } else {
            let r = self.real / m;
            let i = self.imaginary / m;
            (r.mul_add(r, i * i)) * (m * m)
        }
    }

    /// Creates a complex number from polar coordinates.
    ///
    /// # Arguments
    /// - `r`: The abs (radius).
    /// - `theta`: The argument (angle) in radians.
    ///
    /// # Formula
    /// Returns `r * (cos(θ) + i*sin(θ)) = r * e^(iθ)`
    ///
    /// # Examples
    /// ```
    /// use lera::scalar::{Complex, Float};
    ///
    /// let z = Complex::<f64>::from_polar(1.0, f64::FRAC_PI_2);
    /// // z ≈ i (0 + 1i)
    /// assert!(z.real().abs() < 1e-10);
    /// assert!((z.imaginary() - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn from_polar(r: T, theta: T) -> Self {
        let (s, c) = theta.sin_cos();
        Self::new(r * c, r * s)
    }

    /// Converts the complex number's argument from radians to degrees.
    ///
    /// Preserves abs, converts the phase angle to degrees.
    #[inline]
    pub fn to_degrees(self) -> Self {
        let r = self.abs();
        let theta = self.arg().to_degrees();
        Self::from_polar(r, theta.to_radians())
    }

    /// Converts the complex number's argument from degrees to radians.
    ///
    /// Preserves abs, converts the phase angle to radians.
    #[inline]
    pub fn to_radians(self) -> Self {
        let r = self.abs();
        let theta = self.arg().to_radians();
        Self::from_polar(r, theta)
    }

    /// Returns the real part of the complex number.
    #[inline]
    pub const fn real(self) -> T {
        self.real
    }

    /// Returns the imaginary part of the complex number.
    #[inline]
    pub const fn imaginary(self) -> T {
        self.imaginary
    }

    /// Returns true if either component is NaN.
    ///
    /// Uses the IEEE 754 property that NaN != NaN.
    #[inline]
    pub fn is_nan(self) -> bool {
        self.real.is_nan() || self.imaginary.is_nan()
    }

    /// Returns true if either component is infinite.
    #[inline]
    pub fn is_infinite(self) -> bool {
        self.real.is_infinite() || self.imaginary.is_infinite()
    }

    /// Returns true if both components are finite (not NaN and not infinite).
    #[inline]
    pub fn is_finite(self) -> bool {
        self.real.is_finite() && self.imaginary.is_finite()
    }
    /// Returns the additive identity (0 + 0i).
    pub const ZERO: Self = Self::new(T::ZERO, T::ZERO);
    /// Returns the multiplicative identity (1 + 0i).
    pub const ONE: Self = Self::new(T::ONE, T::ZERO);
    /// Returns the constant 2 + 0i.
    pub const TWO: Self = Self::new(T::TWO, T::ZERO);
    /// Returns the imaginary unit (0 + 1i).
    pub const I: Self = Self::new(T::ZERO, T::ONE);
}
