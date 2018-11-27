use std::iter::IntoIterator;

use ::bindings;
use ::{IntegrationResult, Integrator, Real};
use ::ffi::LandingPad;
use ::traits::{IntegrandInput, IntegrandOutput};

use super::{make_gsl_function, GSLIntegrationError, GSLIntegrationWorkspace};

/// Quadrature Adaptive General integration with known singular Points.
/// QAGP applies the same adaptive algorithm as QAGS, with the benefit of known
/// locations of singularities.
///
/// See its documentation
/// [here](https://www.gnu.org/software/gsl/doc/html/integration.html#c.gsl_integration_qagp).
pub struct QAGP {
    singularities: Vec<Real>,
    wkspc: GSLIntegrationWorkspace,
}

fn verify_singular_points<I>(iter: I) -> Option<Vec<Real>>
    where I: IntoIterator<Item=Real> {
    let vec = iter.into_iter().collect::<Vec<Real>>();

    if vec.len() < 2 {
        return None
    }

    for (&a, &b) in vec[..].iter().zip(vec[1..].iter()) {
        if a >= b {
            return None
        }
    }
    Some(vec)
}

impl QAGP {
    /// Creates a new QAGP with enough memory for `nintervals` subintervals,
    /// and integration bounds and singular points defined by `iter`.
    /// The first and last values in `iter` should be the integration bounds,
    /// and all the others locations of known singularities within those
    /// integration bounds.
    /// Returns `None` if there are less than 2 points (because endpoints
    /// would not be defined), or if the singular points are not in ascending
    /// order.
    pub fn new<I>(nintervals: usize, iter: I) -> Option<Self>
        where I: IntoIterator<Item=Real> {
        Some(QAGP {
            singularities: verify_singular_points(iter)?,
            wkspc: GSLIntegrationWorkspace::new(nintervals)
        })
    }

    /// Discards the old workspace and allocates a new one with enough memory
    /// for `nintervals` subintervals.
    pub fn with_nintervals(self, nintervals: usize) -> Self {
        QAGP {
            wkspc: GSLIntegrationWorkspace::new(nintervals),
            ..self
        }
    }

    /// Provides new singular points. As with `QAGP::new()`, returns `None` if
    /// there are less than 2 points, or if the singular points are not in
    /// ascending order.
    pub fn with_points<I>(self, iter: I) -> Option<Self>
        where I: IntoIterator<Item=Real> {
        Some(QAGP {
            singularities: verify_singular_points(iter)?,
            ..self
        })
    }

    pub fn singularities(&self) -> &[Real] {
        &self.singularities[..]
    }
}

impl Integrator for QAGP {
    type Success = IntegrationResult;
    type Failure = GSLIntegrationError;
    fn integrate<A, B, F: FnMut(A) -> B>(&mut self, fun: F, epsrel: Real, epsabs: Real) -> Result<Self::Success, Self::Failure>
        where A: IntegrandInput,
              B: IntegrandOutput
    {
        let (low, high) = (*self.singularities.first().expect("can't be empty"),
                           *self.singularities.last().expect("can't be empty"));
        let singularities = &mut self.singularities[..];
        let mut value: Real = 0.0;
        let mut error: Real = 0.0;

        let mut lp = LandingPad::new(fun);
        let retcode = unsafe {
            let mut gslfn = make_gsl_function(&mut lp, low, high)?;
            bindings::gsl_integration_qagp(&mut gslfn.function,
                                           singularities.as_mut_ptr(),
                                           singularities.len(),
                                           epsabs, epsrel,
                                           self.wkspc.nintervals,
                                           self.wkspc.wkspc,
                                           &mut value,
                                           &mut error)
        };
        lp.maybe_resume_unwind();

        if retcode != bindings::GSL_SUCCESS {
            Err(GSLIntegrationError::GSLError(retcode.into()))
        } else {
            Ok(IntegrationResult {
                value, error
            })
        }
    }
}
