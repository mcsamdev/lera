use std::fmt;

impl fmt::Display for crate::scalar::complex::core::Complex<f64> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.imaginary == 0.0 {
            return write!(f, "{}", self.real);
        }
        write!(f, "({} + {}i)", self.real, self.imaginary)
    }
}
