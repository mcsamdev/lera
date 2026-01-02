#[cfg(test)]
mod tests {
    use crate::scalar::complex::core::Complex;

    const EPSILON: f64 = 1e-10;
    const EPSILON_F32: f32 = 1e-5;

    /// Helper to check if two complex numbers are approximately equal
    fn approx_eq(a: Complex<f64>, b: Complex<f64>, eps: f64) -> bool {
        (a.real - b.real).abs() < eps && (a.imaginary - b.imaginary).abs() < eps
    }

    fn approx_eq_f32(a: Complex<f32>, b: Complex<f32>, eps: f32) -> bool {
        (a.real - b.real).abs() < eps && (a.imaginary - b.imaginary).abs() < eps
    }

    // ==================== Edge Cases Tests ====================

    #[test]
    fn test_div_by_zero() {
        let z = Complex::new(1.0_f64, 1.0);
        let zero = Complex::<f64>::zero();
        let result = z / zero;
        // Division by zero produces NaN (since 0/0 in the denominator calculation)
        // or infinity depending on the values
        assert!(!result.is_finite());
    }

    #[test]
    fn test_ln_zero() {
        let z = Complex::<f64>::zero();
        let result = z.ln();
        // ln(0) should have -∞ real part
        assert_eq!(result.real(), f64::NEG_INFINITY);
    }

    #[test]
    fn test_operations_with_infinity() {
        let inf = Complex::new(f64::INFINITY, 0.0);
        let z = Complex::new(1.0_f64, 1.0);

        // Adding infinity
        let sum = z + inf;
        assert!(sum.is_infinite());

        // Multiplying by infinity
        let prod = z * inf;
        assert!(prod.is_infinite());
    }

    #[test]
    fn test_operations_with_nan() {
        let nan = Complex::new(f64::NAN, 0.0);
        let z = Complex::new(1.0_f64, 1.0);

        // Any operation with NaN should produce NaN
        let sum = z + nan;
        assert!(sum.is_nan());
    }

    #[test]
    fn test_very_small_numbers() {
        let tiny = Complex::new(1e-300_f64, 1e-300);
        let magnitude = tiny.magnitude();
        assert!(magnitude > 0.0);
        assert!(magnitude.is_finite());
    }

    #[test]
    fn test_very_large_numbers() {
        let huge = Complex::new(1e300_f64, 1e300);
        let magnitude = huge.magnitude();
        assert!(magnitude.is_finite());
    }

    // ==================== f32 Tests ====================

    #[test]
    fn test_f32_basic() {
        let z = Complex::new(3.0_f32, 4.0);
        assert!((z.magnitude() - 5.0).abs() < EPSILON_F32);
    }

    #[test]
    fn test_f32_exp_euler() {
        // e^(iπ) = -1
        let z = Complex::new(0.0_f32, std::f32::consts::PI);
        let result = z.exp();
        assert!(approx_eq_f32(result, Complex::new(-1.0, 0.0), EPSILON_F32));
    }

    // ==================== Iterator Tests ====================

    #[test]
    fn test_sum() {
        let numbers = vec![
            Complex::new(1.0_f64, 2.0),
            Complex::new(3.0, 4.0),
            Complex::new(5.0, 6.0),
        ];
        let sum: Complex<f64> = numbers.into_iter().sum();
        assert_eq!(sum, Complex::new(9.0, 12.0));
    }

    #[test]
    fn test_product() {
        let numbers = vec![
            Complex::new(1.0_f64, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(3.0, 0.0),
        ];
        let prod: Complex<f64> = numbers.into_iter().product();
        assert!(approx_eq(prod, Complex::new(6.0, 0.0), EPSILON));
    }

    // ==================== Conversion Tests ====================

    #[test]
    fn test_from_f64() {
        let z: Complex<f64> = 5.0_f64.into();
        assert_eq!(z, Complex::new(5.0, 0.0));
    }

    #[test]
    fn test_from_f32() {
        let z: Complex<f32> = 5.0_f32.into();
        assert_eq!(z, Complex::new(5.0, 0.0));
    }

    #[test]
    fn test_from_i32() {
        let z: Complex<f32> = 5_i32.into();
        assert_eq!(z, Complex::new(5.0, 0.0));
    }

    #[test]
    fn test_from_i64() {
        let z: Complex<f64> = 5_i64.into();
        assert_eq!(z, Complex::new(5.0, 0.0));
    }

    // ==================== Utility Tests ====================

    #[test]
    fn test_is_nan() {
        let z = Complex::new(f64::NAN, 0.0);
        assert!(z.is_nan());

        let z = Complex::new(0.0_f64, f64::NAN);
        assert!(z.is_nan());

        let z = Complex::new(1.0_f64, 1.0);
        assert!(!z.is_nan());
    }

    #[test]
    fn test_is_infinite() {
        let z = Complex::new(f64::INFINITY, 0.0);
        assert!(z.is_infinite());

        let z = Complex::new(0.0_f64, f64::NEG_INFINITY);
        assert!(z.is_infinite());

        let z = Complex::new(1.0_f64, 1.0);
        assert!(!z.is_infinite());
    }

    #[test]
    fn test_is_finite() {
        let z = Complex::new(1.0_f64, 1.0);
        assert!(z.is_finite());

        let z = Complex::new(f64::INFINITY, 0.0);
        assert!(!z.is_finite());

        let z = Complex::new(f64::NAN, 0.0);
        assert!(!z.is_finite());
    }

    #[test]
    fn test_default() {
        let z: Complex<f64> = Complex::default();
        assert_eq!(z, Complex::zero());
    }

    // ==================== Assignment Operator Tests ====================

    #[test]
    fn test_add_assign() {
        let mut z = Complex::new(1.0_f64, 2.0);
        z += Complex::new(3.0, 4.0);
        assert_eq!(z, Complex::new(4.0, 6.0));
    }

    #[test]
    fn test_sub_assign() {
        let mut z = Complex::new(5.0_f64, 7.0);
        z -= Complex::new(3.0, 4.0);
        assert_eq!(z, Complex::new(2.0, 3.0));
    }

    #[test]
    fn test_mul_assign() {
        let mut z = Complex::new(1.0_f64, 1.0);
        z *= Complex::new(1.0, 1.0);
        assert!(approx_eq(z, Complex::new(0.0, 2.0), EPSILON));
    }

    #[test]
    fn test_div_assign() {
        let mut z = Complex::new(4.0_f64, 2.0);
        z /= Complex::new(1.0, 1.0);
        assert!(approx_eq(z, Complex::new(3.0, -1.0), EPSILON));
    }

    #[test]
    fn test_mul_assign_real() {
        let mut z = Complex::new(2.0_f64, 3.0);
        z *= 2.0;
        assert_eq!(z, Complex::new(4.0, 6.0));
    }

    #[test]
    fn test_div_assign_real() {
        let mut z = Complex::new(4.0_f64, 6.0);
        z /= 2.0;
        assert_eq!(z, Complex::new(2.0, 3.0));
    }
}
