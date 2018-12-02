# Integrators
A Rust crate which provides a generic interface for numerical integrators. It includes implementations of this interface for integrators from the GSL and Cuba libraries.


## Examples

This example will integrate a Gaussian over a given range with a GSL integrator. In reality, of course, you should probably find an `erf()` implementation to call instead, but this illustrates its use.

```rust
extern crate integrators;
use integrators::{Integrator, Real};

fn integrate_gaussian(from: f64, to: f64, sigma: f64, mean: f64) -> f64 {
    let normalization = (2f64 * ::std::f64::consts::PI).sqrt() * sigma;
    integrators::gsl::QAG::new(1000)
                      .with_range(from, to)
                      .integrate(|x: Real| {
                          (-((x - mean) / (2f64 * sigma)).powi(2)).exp()
                          / normalization
                          }, 1e-6, 1e-18)
                      .unwrap()
                      .value
}

fn main() {
    let ranges: &'static [(Real, Real)] = &[(-0.3, 0.0), (0.0, 1.5),
                                            (1.5, 3.2), (3.2, 5.9)];
    for &(a, b) in ranges.iter() {
        let integrated = integrate_gaussian(a, b, 1.0, 0.0);
        println!("range: {}, {}", a, b);
        println!("integrated: {}", integrated);
    }
}
```
