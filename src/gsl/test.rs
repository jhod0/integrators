//use std::intrinsics::unchecked_div;
use ::Real;
use ::Integrator;
use super::{GSLIntegrationError, QNG, QAG, QAGS, QAGP, QAWC};

fn nan(_: Real) -> Real {
    ::std::f64::NAN
}

fn infty(_: Real) -> Real {
    ::std::f64::INFINITY
}

fn inv_sq_offset(x: Real, ofs: Real) -> Real {
    1.0 / (x - ofs).powi(2)
}

fn inv_offset(x: Real, ofs: Real) -> Real {
    1.0 / (x - ofs)
}

fn two_args((a, b): (Real, Real)) -> Real {
    a + b
}

fn two_outputs(_: Real) -> (Real, Real) {
    (1.0, 1.0)
}

fn three_inputs_two_outputs((a, b, c): (Real, Real, Real)) -> (Real, Real) {
    (a * b, b * c)
}

fn quadratic_1(x: Real) -> Real {
    x*x - 3f64*x + 17.3
}

fn quadratic_1_integral(range_low: Real, range_high: Real) -> Real {
    let cubic = (range_high.powi(3) - range_low.powi(3)) / 3f64;
    let quadratic = (range_high.powi(2) - range_low.powi(2)) / 2f64;
    cubic - 3f64*quadratic + 17.3*(range_high - range_low)
}

#[test]
fn test_quadratic() {
    let mut qng = QNG::new(0.0, 1.0);
    let mut qag = QAG::new(10);
    let ranges = vec![(0.0, 1.0), (1.0, 2.0), (2.0, 4.0), (-10.0, 30.0)];
    for (low, high) in ranges.into_iter() {
        qng = qng.with_range(low, high);
        qag = qag.with_range(low, high);
        let exp = quadratic_1_integral(low, high);
        let qng_res = qng.integrate(quadratic_1, 1e-3, 1e-6)
                         .expect("should converge");
        let qag_res = qag.integrate(quadratic_1, 1e-3, 1e-6)
                         .expect("should converge");
        assert!((qng_res.value - exp).abs() <= qng_res.error);
        assert!((qag_res.value - exp).abs() <= qag_res.error);
    }
}

#[test]
fn test_invalid_dims() {
    let mut qng = QNG::new(0.0, 1.0);
    assert_eq!(qng.integrate(two_args, 1e-3, 1e-6),
               Err(GSLIntegrationError::InvalidInputDim(2)));
    assert_eq!(qng.integrate(two_outputs, 1e-3, 1e-6),
               Err(GSLIntegrationError::InvalidOutputDim(2)));
    assert_eq!(qng.integrate(three_inputs_two_outputs, 1e-3, 1e-6),
               Err(GSLIntegrationError::InvalidInputDim(3)));
}

#[test]
fn test_error_handling_nan() {
    let mut qng = QNG::new(0.0, 1.0);
    let res = qng.integrate(nan, 1e-6, 1e-10)
                 .expect_err("integration should fail");
    assert_eq!("(GSL) error code 14, description: failed to reach the specified tolerance",
               format!("{}", res));
}

#[test]
fn test_error_handling_pos_inf() {
    let mut qng = QNG::new(0.0, 1.0);
    let res = qng.integrate(infty, 1e-6, 1e-10)
                 .expect_err("integration should fail");
    assert_eq!("(GSL) error code 14, description: failed to reach the specified tolerance",
               format!("{}", res));
}

#[test]
fn test_error_handling_singularity() {
    let mut qag = QAG::new(1000)
                      .with_range(0.0, 1.0);
    let res = qag.integrate(|x| inv_sq_offset(x, 0.5), 1e-6, 1e-10)
                 .expect_err("integration should fail");
    assert_eq!("(GSL) error code 5, description: generic failure",
               format!("{}", res));
}

#[test]
fn test_qags_singularity() {
    let mut qags = QAGS::new(1000);
    let res1 = qags.integrate(|x| inv_offset(x, 0.5), 1e-6, 1e-10)
                   .expect_err("integration should fail");
    println!("{}", res1);
    assert_eq!("(GSL) error code 21, description: singularity or extremely bad function behavior detected",
               format!("{}", res1));

    let res2 = qags.integrate(|x: Real| 1.0 / (1.0 - x).sqrt(), 1e-6, 1e-10)
                   .expect("integration should succeed");
    assert!((res2.value - 2f64).abs() <= res2.error);
}

#[test]
fn test_qagp_singularity() {
    let mut qagp = QAGP::new(1000, vec![0.0, 0.5, 1.0].into_iter()).unwrap();
    let res1 = qagp.integrate(|x| inv_offset(x, 0.5), 1e-6, 1e-10)
                   .expect("integration should work");
    assert!((res1.value - 0f64).abs() <= res1.error);

    let res2 = qagp.with_points(vec![0.0, 1.0].into_iter()).unwrap()
                   .integrate(|x: Real| 1.0 / (1.0 - x).sqrt(), 1e-6, 1e-10)
                   .expect("integration should succeed");
    assert!((res2.value - 2f64).abs() <= res2.error);
}

#[test]
fn test_qawc(){
    let mut qawc = QAWC::new(1000)
        .with_range(-1., 1.)
        .with_singularity(0.0);
    // Trivial example
    let res1 = qawc.integrate(|x: Real| x, 1e-6, 1e-10)
        .expect("integration should succeed");
    assert!((res1.value - 2.0f64).abs() <= res1.error);

    let mut qawc = qawc.with_range(0.0, 2.0)
        .with_singularity(1.0);
    // Analytical value: 2 E SinhIntegral[1] ~ 5.74781168
    let res2 = qawc.integrate(|x: Real| x.exp(), 1e-6, 1e-10)
        .expect("integration should succeed");
    println!("Result: {},   err:{}", res2.value, res2.error);
    assert!((res2.value - 5.747811685312522940687587).abs() <= res2.error);
}
