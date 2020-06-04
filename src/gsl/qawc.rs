use ::bindings;
use ::{IntegrationResult, Integrator, Real};
use ::ffi::LandingPad;
use ::traits::{IntegrandInput, IntegrandOutput};

use super::{make_gsl_function, GSLIntegrationError, GSLIntegrationWorkspace};

/// Quadrature Adaptive Weighed integration for Cauchy principal values
/// Calculates the principal value of a function f(x) with a singular weight
/// 1 / (x - c) centered at a point c
#[derive(Debug, Clone)]
pub struct QAWC {
    range_low: Real,
    range_high: Real,
    c: Real,
    wkspc: GSLIntegrationWorkspace
}

impl QAWC {
    /// Creates a new QAGS with enough memory for `nintervals` subintervals.
    /// This will create a QAG to integrate the range [0, 1].
    /// To change the integration bounds, see `with_range`.
    /// The singularity is set to 0.5 by default and can be changed with `with_singularity`
    pub fn new(nintervals: usize) -> Self{
        QAWC {
            range_low: 0.0,
            range_high: 1.0,
            c: 0.5,
            wkspc: GSLIntegrationWorkspace::new(nintervals)
        }
    }

    /// Discards the old workspace and allocates a new one with enough memory
    /// for `nintervals` subintervals.
    pub fn with_nintervals(self, nintervals: usize) -> Self{
        QAWC {
            wkspc: GSLIntegrationWorkspace::new(nintervals),
            ..self
        }
    }

    pub fn with_range(self, range_low: Real, range_high: Real) -> Self {
        QAWC { range_low, range_high, ..self }
    }

    pub fn with_singularity(self, c: Real) -> Self{
        QAWC {c, ..self}
    }
}

impl Integrator for QAWC {
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
            bindings::gsl_integration_qawc(&mut gslfn.function,
                                           self.range_low, self.range_high,
                                           self.c,
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