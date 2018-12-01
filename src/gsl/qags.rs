use ::bindings;
use ::{IntegrationResult, Integrator, Real};
use ::ffi::LandingPad;
use ::traits::{IntegrandInput, IntegrandOutput};

use super::{make_gsl_function, GSLIntegrationError, GSLIntegrationWorkspace};

/// Quadrature Adaptive General integration with Singularities. Concentrates
/// subintervals around integrable singularities which converge to the solution,
/// using an extrapolation procedure to speed convergence.
#[derive(Debug, Clone)]
pub struct QAGS {
    range_low: Real,
    range_high: Real,
    wkspc: GSLIntegrationWorkspace,
}

impl QAGS {
    /// Creates a new QAGS with enough memory for `nintervals` subintervals.
    /// This will create a QAG to integrate the range [0, 1]. To change the
    /// integration bounds, see `with_range`.
    pub fn new(nintervals: usize) -> Self {
        QAGS {
            range_low: 0.0,
            range_high: 1.0,
            wkspc: GSLIntegrationWorkspace::new(nintervals)
        }
    }

    /// Discards the old workspace and allocates a new one with enough memory
    /// for `nintervals` subintervals.
    pub fn with_nintervals(self, nintervals: usize) -> Self {
        QAGS {
            wkspc: GSLIntegrationWorkspace::new(nintervals),
            ..self
        }
    }

    pub fn with_range(self, range_low: Real, range_high: Real) -> Self {
        QAGS { range_low, range_high, ..self }
    }
}

impl Integrator for QAGS {
    type Success = IntegrationResult;
    type Failure = GSLIntegrationError;
    fn integrate<A, B, F: FnMut(A) -> B>(&mut self, fun: F, epsrel: Real, epsabs: Real) -> Result<Self::Success, Self::Failure>
        where A: IntegrandInput,
              B: IntegrandOutput
    {
        let mut value: Real = 0.0;
        let mut error: Real = 0.0;

        let mut lp = LandingPad::new(fun);
        let retcode = unsafe {
            let mut gslfn = make_gsl_function(&mut lp, self.range_low, self.range_high)?;
            bindings::gsl_integration_qags(&mut gslfn.function,
                                           self.range_low, self.range_high,
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
