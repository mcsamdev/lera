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

    // ==================== Construction Tests ====================

    #[test]
    fn test_new() {
        let z = Complex::new(3.0_f64, 4.0);
        assert_eq!(z.real(), 3.0);
        assert_eq!(z.imaginary(), 4.0);
    }

    #[test]
    fn test_zero() {
        let z = Complex::<f64>::zero();
        assert_eq!(z.real(), 0.0);
        assert_eq!(z.imaginary(), 0.0);
    }

    #[test]
    fn test_one() {
        let z = Complex::<f64>::one();
        assert_eq!(z.real(), 1.0);
        assert_eq!(z.imaginary(), 0.0);
    }

    #[test]
    fn test_i() {
        let z = Complex::<f64>::i();
        assert_eq!(z.real(), 0.0);
        assert_eq!(z.imaginary(), 1.0);
    }

    #[test]
    fn test_from_polar() {
        // 1 at angle 0 should be (1, 0)
        let z = Complex::<f64>::from_polar(1.0, 0.0);
        assert!(approx_eq(z, Complex::new(1.0, 0.0), EPSILON));

        // 1 at angle π/2 should be (0, 1)
        let z = Complex::<f64>::from_polar(1.0, std::f64::consts::FRAC_PI_2);
        assert!(approx_eq(z, Complex::new(0.0, 1.0), EPSILON));

        // 1 at angle π should be (-1, 0)
        let z = Complex::<f64>::from_polar(1.0, std::f64::consts::PI);
        assert!(approx_eq(z, Complex::new(-1.0, 0.0), EPSILON));

        // √2 at angle π/4 should be (1, 1)
        let z = Complex::<f64>::from_polar(std::f64::consts::SQRT_2, std::f64::consts::FRAC_PI_4);
        assert!(approx_eq(z, Complex::new(1.0, 1.0), EPSILON));
    }

    // ==================== Basic Properties Tests ====================

    #[test]
    fn test_conjugate() {
        let z = Complex::new(3.0_f64, 4.0);
        let conj = z.conjugate();
        assert_eq!(conj.real(), 3.0);
        assert_eq!(conj.imaginary(), -4.0);

        // Conjugate of conjugate is original
        assert_eq!(conj.conjugate(), z);
    }

    #[test]
    fn test_magnitude() {
        // 3-4-5 triangle
        let z = Complex::new(3.0_f64, 4.0);
        assert!((z.magnitude() - 5.0).abs() < EPSILON);

        // Pure real
        let z = Complex::new(5.0_f64, 0.0);
        assert!((z.magnitude() - 5.0).abs() < EPSILON);

        // Pure imaginary
        let z = Complex::new(0.0_f64, 5.0);
        assert!((z.magnitude() - 5.0).abs() < EPSILON);

        // Zero
        let z = Complex::<f64>::zero();
        assert_eq!(z.magnitude(), 0.0);
    }

    #[test]
    fn test_norm_sqr() {
        let z = Complex::new(3.0_f64, 4.0);
        assert!((z.norm_sqr() - 25.0).abs() < EPSILON);
    }

    #[test]
    fn test_arg() {
        // Positive real axis
        let z = Complex::new(1.0_f64, 0.0);
        assert!(z.arg().abs() < EPSILON);

        // Positive imaginary axis
        let z = Complex::new(0.0_f64, 1.0);
        assert!((z.arg() - std::f64::consts::FRAC_PI_2).abs() < EPSILON);

        // Negative real axis
        let z = Complex::new(-1.0_f64, 0.0);
        assert!((z.arg() - std::f64::consts::PI).abs() < EPSILON);

        // First quadrant (45°)
        let z = Complex::new(1.0_f64, 1.0);
        assert!((z.arg() - std::f64::consts::FRAC_PI_4).abs() < EPSILON);
    }

    // ==================== Arithmetic Tests ====================

    #[test]
    fn test_add() {
        let a = Complex::new(1.0_f64, 2.0);
        let b = Complex::new(3.0, 4.0);
        let sum = a + b;
        assert_eq!(sum, Complex::new(4.0, 6.0));
    }

    #[test]
    fn test_add_real() {
        let z = Complex::new(1.0_f64, 2.0);
        let sum = z + 5.0;
        assert_eq!(sum, Complex::new(6.0, 2.0));
    }

    #[test]
    fn test_sub() {
        let a = Complex::new(5.0_f64, 7.0);
        let b = Complex::new(3.0, 4.0);
        let diff = a - b;
        assert_eq!(diff, Complex::new(2.0, 3.0));
    }

    #[test]
    fn test_mul() {
        // (1 + 2i)(3 + 4i) = 3 + 4i + 6i + 8i² = 3 + 10i - 8 = -5 + 10i
        let a = Complex::new(1.0_f64, 2.0);
        let b = Complex::new(3.0, 4.0);
        let prod = a * b;
        assert!(approx_eq(prod, Complex::new(-5.0, 10.0), EPSILON));
    }

    #[test]
    fn test_mul_i() {
        // (3 + 4i) * i = -4 + 3i
        let z = Complex::new(3.0_f64, 4.0);
        let result = z.mul_i();
        assert!(approx_eq(result, Complex::new(-4.0, 3.0), EPSILON));
    }

    #[test]
    fn test_mul_i_neg() {
        // (3 + 4i) * (-i) = 4 - 3i
        let z = Complex::new(3.0_f64, 4.0);
        let result = z.mul_i_neg();
        assert!(approx_eq(result, Complex::new(4.0, -3.0), EPSILON));
    }

    #[test]
    fn test_div() {
        // (1 + 2i) / (3 + 4i)
        let a = Complex::new(1.0_f64, 2.0);
        let b = Complex::new(3.0, 4.0);
        let quot = a / b;
        // Expected: (1+2i)(3-4i) / 25 = (3 - 4i + 6i + 8) / 25 = (11 + 2i) / 25
        let expected = Complex::new(11.0 / 25.0, 2.0 / 25.0);
        assert!(approx_eq(quot, expected, EPSILON));

        // Verify: (a / b) * b ≈ a
        let check = quot * b;
        assert!(approx_eq(check, a, EPSILON));
    }

    #[test]
    fn test_recip() {
        let z = Complex::new(3.0_f64, 4.0);
        let inv = z.recip();

        // z * (1/z) should equal 1
        let product = z * inv;
        assert!(approx_eq(product, Complex::one(), EPSILON));
    }

    #[test]
    fn test_neg() {
        let z = Complex::new(3.0_f64, 4.0);
        let neg = -z;
        assert_eq!(neg, Complex::new(-3.0, -4.0));
    }

    // ==================== Power and Root Tests ====================

    #[test]
    fn test_powi_positive() {
        let z = Complex::new(1.0_f64, 1.0);

        // z^0 = 1
        assert!(approx_eq(z.powi(0), Complex::one(), EPSILON));

        // z^1 = z
        assert!(approx_eq(z.powi(1), z, EPSILON));

        // z^2 = (1+i)² = 2i
        let z2 = z.powi(2);
        assert!(approx_eq(z2, Complex::new(0.0, 2.0), EPSILON));

        // z^4 = (2i)² = -4
        let z4 = z.powi(4);
        assert!(approx_eq(z4, Complex::new(-4.0, 0.0), EPSILON));
    }

    #[test]
    fn test_powi_negative() {
        let z = Complex::new(1.0_f64, 1.0);

        // z^(-1) = 1/z
        let z_neg1 = z.powi(-1);
        let expected = z.recip();
        assert!(approx_eq(z_neg1, expected, EPSILON));

        // z^(-2) should be (1/z)²
        let z_neg2 = z.powi(-2);
        let expected = z.recip().powi(2);
        assert!(approx_eq(z_neg2, expected, EPSILON));
    }

    #[test]
    fn test_powi_i_powers() {
        let i = Complex::<f64>::i();

        // i^1 = i
        assert!(approx_eq(i.powi(1), Complex::new(0.0, 1.0), EPSILON));

        // i^2 = -1
        assert!(approx_eq(i.powi(2), Complex::new(-1.0, 0.0), EPSILON));

        // i^3 = -i
        assert!(approx_eq(i.powi(3), Complex::new(0.0, -1.0), EPSILON));

        // i^4 = 1
        assert!(approx_eq(i.powi(4), Complex::new(1.0, 0.0), EPSILON));
    }

    #[test]
    fn test_sqrt() {
        // sqrt(1) = 1
        let z = Complex::new(1.0_f64, 0.0);
        let root = z.sqrt();
        assert!(approx_eq(root, Complex::new(1.0, 0.0), EPSILON));

        // sqrt(-1) = i
        let z = Complex::new(-1.0_f64, 0.0);
        let root = z.sqrt();
        assert!(approx_eq(root, Complex::new(0.0, 1.0), EPSILON));

        // sqrt(i) = (1+i)/√2
        let z = Complex::<f64>::i();
        let root = z.sqrt();
        let expected = Complex::new(
            1.0 / std::f64::consts::SQRT_2,
            1.0 / std::f64::consts::SQRT_2,
        );
        assert!(approx_eq(root, expected, EPSILON));

        // Verify: sqrt(z)² = z
        let z = Complex::new(3.0_f64, 4.0);
        let root = z.sqrt();
        let check = root * root;
        assert!(approx_eq(check, z, EPSILON));
    }

    #[test]
    fn test_sqrt_zero() {
        let z = Complex::<f64>::zero();
        let root = z.sqrt();
        assert!(approx_eq(root, Complex::zero(), EPSILON));
    }

    #[test]
    fn test_cbrt() {
        // cbrt(1) = 1
        let z = Complex::new(1.0_f64, 0.0);
        let root = z.cbrt();
        assert!(approx_eq(root, Complex::new(1.0, 0.0), EPSILON));

        // cbrt(8) = 2
        let z = Complex::new(8.0_f64, 0.0);
        let root = z.cbrt();
        assert!(approx_eq(root, Complex::new(2.0, 0.0), EPSILON));

        // Verify: cbrt(z)³ ≈ z
        let z = Complex::new(3.0_f64, 4.0);
        let root = z.cbrt();
        let check = root.powi(3);
        assert!(approx_eq(check, z, 1e-8));
    }

    #[test]
    fn test_powf() {
        // 2^0.5 = √2
        let z = Complex::new(2.0_f64, 0.0);
        let result = z.powf(0.5);
        assert!((result.real() - std::f64::consts::SQRT_2).abs() < EPSILON);
        assert!(result.imaginary().abs() < EPSILON);

        // Verify z^n using powf matches powi for integer n
        let z = Complex::new(2.0_f64, 1.0);
        let pow_i = z.powi(3);
        let pow_f = z.powf(3.0);
        assert!(approx_eq(pow_i, pow_f, 1e-8));
    }

    // ==================== Exponential and Logarithm Tests ====================

    #[test]
    fn test_exp() {
        // e^0 = 1
        let z = Complex::<f64>::zero();
        let result = z.exp();
        assert!(approx_eq(result, Complex::one(), EPSILON));

        // e^1 = e
        let z = Complex::new(1.0_f64, 0.0);
        let result = z.exp();
        assert!((result.real() - std::f64::consts::E).abs() < EPSILON);
        assert!(result.imaginary().abs() < EPSILON);

        // e^(iπ) = -1 (Euler's identity)
        let z = Complex::new(0.0_f64, std::f64::consts::PI);
        let result = z.exp();
        assert!(approx_eq(result, Complex::new(-1.0, 0.0), EPSILON));

        // e^(iπ/2) = i
        let z = Complex::new(0.0_f64, std::f64::consts::FRAC_PI_2);
        let result = z.exp();
        assert!(approx_eq(result, Complex::new(0.0, 1.0), EPSILON));
    }

    #[test]
    fn test_exp2() {
        // 2^0 = 1
        let z = Complex::<f64>::zero();
        let result = z.exp2();
        assert!(approx_eq(result, Complex::one(), EPSILON));

        // 2^1 = 2
        let z = Complex::new(1.0_f64, 0.0);
        let result = z.exp2();
        assert!((result.real() - 2.0).abs() < EPSILON);

        // 2^3 = 8
        let z = Complex::new(3.0_f64, 0.0);
        let result = z.exp2();
        assert!((result.real() - 8.0).abs() < EPSILON);
    }

    #[test]
    fn test_ln() {
        // ln(1) = 0
        let z = Complex::new(1.0_f64, 0.0);
        let result = z.ln();
        assert!(approx_eq(result, Complex::zero(), EPSILON));

        // ln(e) = 1
        let z = Complex::new(std::f64::consts::E, 0.0);
        let result = z.ln();
        assert!(approx_eq(result, Complex::new(1.0, 0.0), EPSILON));

        // ln(-1) = iπ
        let z = Complex::new(-1.0_f64, 0.0);
        let result = z.ln();
        assert!(approx_eq(
            result,
            Complex::new(0.0, std::f64::consts::PI),
            EPSILON
        ));

        // ln(i) = iπ/2
        let z = Complex::<f64>::i();
        let result = z.ln();
        assert!(approx_eq(
            result,
            Complex::new(0.0, std::f64::consts::FRAC_PI_2),
            EPSILON
        ));

        // exp(ln(z)) = z
        let z = Complex::new(3.0_f64, 4.0);
        let result = z.ln().exp();
        assert!(approx_eq(result, z, EPSILON));
    }

    #[test]
    fn test_log2() {
        // log2(1) = 0
        let z = Complex::new(1.0_f64, 0.0);
        let result = z.log2();
        assert!(approx_eq(result, Complex::zero(), EPSILON));

        // log2(2) = 1
        let z = Complex::new(2.0_f64, 0.0);
        let result = z.log2();
        assert!(approx_eq(result, Complex::new(1.0, 0.0), EPSILON));

        // log2(8) = 3
        let z = Complex::new(8.0_f64, 0.0);
        let result = z.log2();
        assert!(approx_eq(result, Complex::new(3.0, 0.0), EPSILON));
    }

    #[test]
    fn test_log10() {
        // log10(1) = 0
        let z = Complex::new(1.0_f64, 0.0);
        let result = z.log10();
        assert!(approx_eq(result, Complex::zero(), EPSILON));

        // log10(10) = 1
        let z = Complex::new(10.0_f64, 0.0);
        let result = z.log10();
        assert!(approx_eq(result, Complex::new(1.0, 0.0), EPSILON));

        // log10(100) = 2
        let z = Complex::new(100.0_f64, 0.0);
        let result = z.log10();
        assert!(approx_eq(result, Complex::new(2.0, 0.0), EPSILON));
    }

    #[test]
    fn test_log_arbitrary_base() {
        // log_3(27) = 3
        let z = Complex::new(27.0_f64, 0.0);
        let result = z.log(3.0);
        assert!(approx_eq(result, Complex::new(3.0, 0.0), EPSILON));
    }

    // ==================== Trigonometric Tests ====================

    #[test]
    fn test_sin() {
        // sin(0) = 0
        let z = Complex::<f64>::zero();
        let result = z.sin();
        assert!(approx_eq(result, Complex::zero(), EPSILON));

        // sin(π/2) = 1
        let z = Complex::new(std::f64::consts::FRAC_PI_2, 0.0);
        let result = z.sin();
        assert!(approx_eq(result, Complex::new(1.0, 0.0), EPSILON));

        // sin(π) ≈ 0
        let z = Complex::new(std::f64::consts::PI, 0.0);
        let result = z.sin();
        assert!(approx_eq(result, Complex::zero(), EPSILON));
    }

    #[test]
    fn test_cos() {
        // cos(0) = 1
        let z = Complex::<f64>::zero();
        let result = z.cos();
        assert!(approx_eq(result, Complex::one(), EPSILON));

        // cos(π/2) ≈ 0
        let z = Complex::new(std::f64::consts::FRAC_PI_2, 0.0);
        let result = z.cos();
        assert!(approx_eq(result, Complex::zero(), EPSILON));

        // cos(π) = -1
        let z = Complex::new(std::f64::consts::PI, 0.0);
        let result = z.cos();
        assert!(approx_eq(result, Complex::new(-1.0, 0.0), EPSILON));
    }

    #[test]
    fn test_tan() {
        // tan(0) = 0
        let z = Complex::<f64>::zero();
        let result = z.tan();
        assert!(approx_eq(result, Complex::zero(), EPSILON));

        // tan(π/4) = 1
        let z = Complex::new(std::f64::consts::FRAC_PI_4, 0.0);
        let result = z.tan();
        assert!(approx_eq(result, Complex::new(1.0, 0.0), EPSILON));
    }

    #[test]
    fn test_sin_cos_identity() {
        // sin²(z) + cos²(z) = 1
        let z = Complex::new(1.0_f64, 2.0);
        let sin_z = z.sin();
        let cos_z = z.cos();
        let sum = sin_z * sin_z + cos_z * cos_z;
        assert!(approx_eq(sum, Complex::one(), EPSILON));
    }

    #[test]
    fn test_asin() {
        // asin(0) = 0
        let z = Complex::<f64>::zero();
        let result = z.asin();
        assert!(approx_eq(result, Complex::zero(), EPSILON));

        // asin(1) = π/2
        let z = Complex::new(1.0_f64, 0.0);
        let result = z.asin();
        assert!(approx_eq(
            result,
            Complex::new(std::f64::consts::FRAC_PI_2, 0.0),
            EPSILON
        ));

        // sin(asin(z)) = z
        let z = Complex::new(0.5_f64, 0.3);
        let result = z.asin().sin();
        assert!(approx_eq(result, z, 1e-8));
    }

    #[test]
    fn test_acos() {
        // acos(1) = 0
        let z = Complex::new(1.0_f64, 0.0);
        let result = z.acos();
        assert!(approx_eq(result, Complex::zero(), EPSILON));

        // cos(acos(z)) = z
        let z = Complex::new(0.5_f64, 0.3);
        let result = z.acos().cos();
        assert!(approx_eq(result, z, 1e-8));
    }

    #[test]
    fn test_atan() {
        // atan(0) = 0
        let z = Complex::<f64>::zero();
        let result = z.atan();
        assert!(approx_eq(result, Complex::zero(), EPSILON));

        // atan(1) = π/4
        let z = Complex::new(1.0_f64, 0.0);
        let result = z.atan();
        assert!(approx_eq(
            result,
            Complex::new(std::f64::consts::FRAC_PI_4, 0.0),
            EPSILON
        ));

        // tan(atan(z)) = z
        let z = Complex::new(0.5_f64, 0.3);
        let result = z.atan().tan();
        assert!(approx_eq(result, z, 1e-8));
    }

    // ==================== Hyperbolic Tests ====================

    #[test]
    fn test_sinh() {
        // sinh(0) = 0
        let z = Complex::<f64>::zero();
        let result = z.sinh();
        assert!(approx_eq(result, Complex::zero(), EPSILON));

        // sinh(1) for real
        let z = Complex::new(1.0_f64, 0.0);
        let result = z.sinh();
        assert!((result.real() - 1.0_f64.sinh()).abs() < EPSILON);
    }

    #[test]
    fn test_cosh() {
        // cosh(0) = 1
        let z = Complex::<f64>::zero();
        let result = z.cosh();
        assert!(approx_eq(result, Complex::one(), EPSILON));

        // cosh(1) for real
        let z = Complex::new(1.0_f64, 0.0);
        let result = z.cosh();
        assert!((result.real() - 1.0_f64.cosh()).abs() < EPSILON);
    }

    #[test]
    fn test_tanh() {
        // tanh(0) = 0
        let z = Complex::<f64>::zero();
        let result = z.tanh();
        assert!(approx_eq(result, Complex::zero(), EPSILON));
    }

    #[test]
    fn test_cosh_sinh_identity() {
        // cosh²(z) - sinh²(z) = 1
        let z = Complex::new(1.0_f64, 2.0);
        let sinh_z = z.sinh();
        let cosh_z = z.cosh();
        let diff = cosh_z * cosh_z - sinh_z * sinh_z;
        assert!(approx_eq(diff, Complex::one(), EPSILON));
    }

    #[test]
    fn test_asinh() {
        // asinh(0) = 0
        let z = Complex::<f64>::zero();
        let result = z.asinh();
        assert!(approx_eq(result, Complex::zero(), EPSILON));

        // sinh(asinh(z)) = z
        let z = Complex::new(0.5_f64, 0.3);
        let result = z.asinh().sinh();
        assert!(approx_eq(result, z, 1e-8));
    }

    #[test]
    fn test_acosh() {
        // acosh(1) = 0
        let z = Complex::new(1.0_f64, 0.0);
        let result = z.acosh();
        assert!(approx_eq(result, Complex::zero(), EPSILON));

        // cosh(acosh(z)) = z
        let z = Complex::new(2.0_f64, 0.3);
        let result = z.acosh().cosh();
        assert!(approx_eq(result, z, 1e-8));
    }

    #[test]
    fn test_atanh() {
        // atanh(0) = 0
        let z = Complex::<f64>::zero();
        let result = z.atanh();
        assert!(approx_eq(result, Complex::zero(), EPSILON));

        // tanh(atanh(z)) = z for |z| < 1
        // Use a value well within the domain for numerical stability
        let z = Complex::new(0.3_f64, 0.2);
        let result = z.atanh().tanh();
        assert!(approx_eq(result, z, 1e-10));
    }
}
