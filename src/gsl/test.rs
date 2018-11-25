use ::Real;
use ::Integrator;
use super::{GSLIntegrationError, QNG, QAG};

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
