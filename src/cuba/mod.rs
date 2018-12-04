//! Wrappers for Cuba integration routines. See Cuba documentation
//! [here](http://www.feynarts.de/cuba/) for details and installation
//! instructions.
//!
//! This module, and all its re-exports from the top-level `integrators`
//! module, are gated with the `cuba` feature. So, if you don't want to use
//! these wrappers, or don't have Cuba installed, just turn off that feature.
//!
//!
//! # Example
//!
//! This example computes the volume of a sphere by direct integration over 3 dimensions, and
//! compares to the expected result of 4/3 pi r^3.
//!
//! ```
//! use std::f64::consts::PI;
//! use integrators::{Integrator, Real, Real3};
//! use integrators::cuba::{Cuhre, IntegrationRange};
//!
//! // An integral over a sphere is:
//! // \int f(\theta, \phi, r) r^2 sin(\theta) d\theta d\phi dr
//! // for r from 0 to the sphere's radius, \theta from 0 to pi,
//! // and \phi from 0 to 2pi
//!
//! fn volume(radius: Real) -> Real {
//!     let mut cuhre = Cuhre::new(999999);
//!     let r_range = IntegrationRange::new(0.0, radius);
//!     let theta_range = IntegrationRange::new(0.0, PI);
//!     let phi_range = IntegrationRange::new(0.0, 2.0*PI);
//!
//!     let results = cuhre.integrate(|(rs, thetas, phis): Real3| {
//!                 let r = r_range.transform(rs);
//!                 let theta = theta_range.transform(thetas);
//!                 let phi = phi_range.transform(phis);
//!                 r * r * theta.sin()
//!                 * r_range.jacobian()
//!                 * theta_range.jacobian()
//!                 * phi_range.jacobian()
//!             }, 1e-5, 1e-18).unwrap();
//!     results.results[0].value
//! }
//!
//! fn exact_volume(radius: Real) -> Real {
//!     4.0/3.0 * PI * radius.powi(3)
//! }
//!
//! for &r in [0.3, 1.0, 5.0, 10.0].into_iter() {
//!     let ex = exact_volume(r);
//!     let calc = volume(r);
//!     assert!((calc - ex).abs() < ex*1e-5);
//! }
//! ```

use std::{error, fmt, slice, vec};
use std::convert::From;
use std::os::raw::{c_int, c_longlong, c_void};

use super::traits::{IntegrandInput, IntegrandOutput};
use super::{IntegrationResult, Real};
use super::ffi::LandingPad;

mod cuhre;
pub use self::cuhre::Cuhre;

mod suave;
pub use self::suave::Suave;

mod vegas;
pub use self::vegas::Vegas;

unsafe extern "C"
fn cuba_integrand<A, B, F>(ndim: *const c_int,
                           x: *const Real,
                           ncomp: *const c_int,
                           f: *mut Real,
                           userdata: *mut c_void) -> c_int
    where A: IntegrandInput,
          B: IntegrandOutput,
          F: FnMut(A) -> B
{
    let fnptr = userdata as *mut LandingPad<A, B, F>;
    let lp: &mut LandingPad<A, B, F> = &mut *fnptr;

    let args = slice::from_raw_parts(x, *ndim as usize);
    let output = slice::from_raw_parts_mut(f, *ncomp as usize);

    match lp.try_call(args, output) {
        Ok(_) => 0,
        // -999 is special `abort` code to Cuba
        Err(_) => -999,
    }
}

/// Since Cuba integrates on the unit hypercube, it is convenient to have a
/// helper to convert into a different integration range.
#[derive(Debug, PartialEq)]
pub struct IntegrationRange {
    start: Real,
    length: Real,
}

impl IntegrationRange {
    /// Create a new range from start to end.
    pub fn new(start: Real, end: Real) -> Self {
        IntegrationRange {
            start,
            length: end - start
        }
    }

    /// Scale Cuba's dimension to this integration range. I.e., convert from
    /// [0,1] to [`start`, `end`].
    ///
    /// # Panics
    /// If x is not between 0 and 1 (inclusive), panics.
    pub fn transform(&self, x: Real) -> Real {
        assert!((x >= 0.0) & (x <= 1.0));
        self.start + x * self.length
    }

    /// Scale factor.
    pub fn jacobian(&self) -> Real {
        self.length
    }
}

/// The random number generator source for Cuba's Monte Carlo algorithms. Refer
/// to Cuba's docs for details and pros/cons of each.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum RandomNumberSource {
    Sobol,
    MersenneTwister,
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct CubaIntegrationResult {
    /// The integration result.
    pub value: Real,
    /// An error estimate on the integration result.
    pub error: Real,
    /// How much trust Cuba puts in the error estimate. Will be a number
    /// between 0 and 1, where 0 is good, and 1 is bad.
    pub prob: Real
}

#[derive(Clone, Debug, PartialEq)]
pub struct CubaIntegrationResults {
    /// The number of subintervals used in the integration. Vegas does not
    /// provide this, and so returns `None` for this field.
    pub nregions: Option<c_int>,
    /// The number of evaluations used.
    pub neval: c_longlong,
    /// Integration results, a vector of the same length as the integrand's
    /// output dimensions.
    pub results: Vec<CubaIntegrationResult>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum CubaError {
    /// The integrand input's dimensions are not supported by the given
    /// algorithm. The name of the algorithm and the number of dimensions
    /// attempted are given.
    BadDim(&'static str, usize),
    /// The integrand output's dimensions are not supported by the given
    /// algorithm. The name of the algorithm and the number of dimensions
    /// attempted are given.
    BadComp(&'static str, usize),
    /// The integration did not converge. Though the results did not reach
    /// the desired uncertainty, they still might be useful, and so are
    /// provided.
    DidNotConverge(CubaIntegrationResults),
}

impl fmt::Display for CubaError {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        use self::CubaError::*;
        match &self {
            &BadDim(name, ndim) => {
                write!(fmt, "invalid number of dimensions for algorithm {}: {}",
                       name, ndim)
            },
            &BadComp(name, ncomp) => {
                write!(fmt, "invalid number of outputs for algorithm {}: {}",
                       name, ncomp)
            },
            &DidNotConverge(_) => write!(fmt, "integral did not converge")
        }
    }
}

impl error::Error for CubaError {}

pub struct CubaResultsIter {
    iter: vec::IntoIter<CubaIntegrationResult>
}

impl Iterator for CubaResultsIter {
    type Item = IntegrationResult;
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|CubaIntegrationResult { value, error, .. }| {
            IntegrationResult {
                value, error
            }
        })
    }
}

impl From<Vec<CubaIntegrationResult>> for CubaResultsIter {
    fn from(v: Vec<CubaIntegrationResult>) -> Self {
        CubaResultsIter {
            iter: v.into_iter()
        }
    }
}

impl super::traits::IntegrationResults for CubaIntegrationResults {
    type Iterator = CubaResultsIter;
    fn results(self) -> CubaResultsIter {
        From::from(self.results)
    }
}
