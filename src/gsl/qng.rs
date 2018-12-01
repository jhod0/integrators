use ::bindings;
use ::ffi::LandingPad;
use ::traits::{IntegrandInput, IntegrandOutput};
use ::{IntegrationResult, Integrator, Real};
use super::{make_gsl_function, GSLIntegrationError};

/// Quadrature Non-adaptive General-use integrator. Iteratively
/// applies Gauss-Kronrod quadrature rules of successively higher
/// order until a desired precision is reached. If such precision
/// is not reached by the 87-point quadrature rule, the integration
/// fails.
///
/// See GSL Docs [here](https://www.gnu.org/software/gsl/doc/html/integration.html#c.gsl_integration_qng).
#[derive(Debug, Copy, Clone)]
pub struct QNG {
    range_low: Real,
    range_high: Real
}

impl QNG {
    /// Creates a new QNG integrator which will integrate a 1-dimensional
    /// function from `range_low` to `range_high`.
    pub fn new(range_low: Real, range_high: Real) -> Self {
        QNG { range_low, range_high }
    }

    /// Builder pattern to change the integration range of the QNG.
    pub fn with_range(self, range_low: Real, range_high: Real) -> Self {
        QNG { range_low, range_high }
    }
}

impl Integrator for QNG {
    type Success = IntegrationResult;
    type Failure = GSLIntegrationError;
    fn integrate<A, B, F: FnMut(A) -> B>(&mut self, fun: F, epsrel: Real, epsabs: Real) -> Result<Self::Success, Self::Failure>
        where A: IntegrandInput,
              B: IntegrandOutput
    {
        let mut value: Real = 0.0;
        let mut error: Real = 0.0;
        let mut neval: usize = 0;

        let mut lp = LandingPad::new(fun);
        let retcode = unsafe {
            let mut gslfn = make_gsl_function(&mut lp, self.range_low, self.range_high)?;
            bindings::gsl_integration_qng(&mut gslfn.function,
                                          self.range_low, self.range_high,
                                          epsabs, epsrel,
                                          &mut value,
                                          &mut error,
                                          &mut neval)
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
