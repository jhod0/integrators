# Integrators
A Rust crate which provides a generic interface for numerical integrators. It includes implementations of this interface for integrators from the GSL and Cuba libraries.


## GSL Wrappers

On Ubuntu-based systems, these wrappers should work with `libgsl-dev`. If you do not wish to use these wrappers, you can disable them by disabling the `gsl` feature.
These algorithms can only integrate one dimensional integrals.

See the GSL docs [here](https://www.gnu.org/software/gsl/doc/html/integration.html#) for a complete list of integration algorithms. Currently, only four wrappers
are implemented:

1. [QAG](https://www.gnu.org/software/gsl/doc/html/integration.html#qag-adaptive-integration) is a general, adaptive integration algorithm which should work for most well-behaved functions.
1. [QNG](https://www.gnu.org/software/gsl/doc/html/integration.html#qng-non-adaptive-gauss-kronrod-integration) is a similarly general, *non*-adaptive algorithm, which applies a series of fixed-order quadrature rules. This algorithm requires less overhead than QAG, and so may provide a performance boost for functions which are known to be easily integrable.
1. [QAGS](https://www.gnu.org/software/gsl/doc/html/integration.html#qags-adaptive-integration-with-singularities) is an adaptive, general algorithm which can handle some kinds of singularities.
1. [QAGP](https://www.gnu.org/software/gsl/doc/html/integration.html#qagp-adaptive-integration-with-known-singular-points) is the same algorithm as QAGS, but requires the user to provide a list of known locations of singularities.

I will add wrappers for more functions as I go.

## Cuba Wrappers

Cuba is a suite of advanced multidimensional numerical integration algorithms, including both Monte Carlo and deterministic
methods. See its documentation [here](http://www.feynarts.de/cuba) for details and installation instructions.
If you do not wish to use these wrappers, you can disable them by disabling the `cuba` feature.

Cuba has four algorithms: Vegas, Suave, Cuhre, and Divonne. Currently, wrappers are only implemented for the first three; the last
has a number of extra arguments for finding singularities which I have not addressed yet.

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
