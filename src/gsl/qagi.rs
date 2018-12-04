use std::convert::Into;

use ::bindings;
use ::{IntegrationResult, Integrator, Real};
use ::ffi::LandingPad;
use ::traits::{IntegrandInput, IntegrandOutput};

use super::{make_gsl_function, GSLIntegrationError, GSLIntegrationWorkspace};

/// Quadrature Adaptive General integration for Infinite intervals. It applies
/// the QAGS algorithm to a transformation of the input integral, such that
/// integrating the transformed function from 0 to 1 yields the infinite
/// integral of the original function.
///
/// See GSL docs
/// [here](https://www.gnu.org/software/gsl/doc/html/integration.html#c.gsl_integration_qagi).
#[derive(Debug, Clone)]
pub struct QAGI {
    wkspc: GSLIntegrationWorkspace,
}

/// Quadrature Adaptive General integration for Infinite intervals, over the
/// semi-infinite interval `(a, +inf)`.
#[derive(Debug, Clone)]
pub struct QAGIU {
    wkspc: GSLIntegrationWorkspace,
    lower_bound: Real,
}

/// Quadrature Adaptive General integration for Infinite intervals, over the
/// semi-infinite interval `(-inf, b)`.
#[derive(Debug, Clone)]
pub struct QAGIL {
    wkspc: GSLIntegrationWorkspace,
    upper_bound: Real,
}

impl QAGI {
    /// Creates a new QAGI with enough memory for `nintervals` subintervals.
    pub fn new(nintervals: usize) -> Self {
        QAGI {
            wkspc: GSLIntegrationWorkspace::new(nintervals)
        }
    }

    /// Discards the old workspace and allocates a new one with enough memory
    /// for `nintervals` subintervals.
    pub fn with_nintervals(self, nintervals: usize) -> Self {
        QAGI {
            wkspc: GSLIntegrationWorkspace::new(nintervals),
            ..self
        }
    }
}

impl QAGIU {
    /// Creates a new QAGIU with enough memory for `nintervals` subintervals,
    /// which will integrate from `lower_bound` to +infinity.
    pub fn new(nintervals: usize, lower_bound: Real) -> Self {
        QAGIU {
            wkspc: GSLIntegrationWorkspace::new(nintervals),
            lower_bound
        }
    }

    /// Discards the old workspace and allocates a new one with enough memory
    /// for `nintervals` subintervals.
    pub fn with_nintervals(self, nintervals: usize) -> Self {
        QAGIU {
            wkspc: GSLIntegrationWorkspace::new(nintervals),
            ..self
        }
    }
}

impl QAGIL {
    /// Creates a new QAGIL with enough memory for `nintervals` subintervals,
    /// which will integrate from -infinity to `upper_bound`.
    pub fn new(nintervals: usize, upper_bound: Real) -> Self {
        QAGIL {
            wkspc: GSLIntegrationWorkspace::new(nintervals),
            upper_bound
        }
    }

    /// Discards the old workspace and allocates a new one with enough memory
    /// for `nintervals` subintervals.
    pub fn with_nintervals(self, nintervals: usize) -> Self {
        QAGIL {
            wkspc: GSLIntegrationWorkspace::new(nintervals),
            ..self
        }
    }
}

impl Integrator for QAGI {
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
            let mut gslfn = make_gsl_function(&mut lp, -1.0, 1.0)?;
            bindings::gsl_integration_qagi(&mut gslfn.function,
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

impl Integrator for QAGIU {
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
            let mut gslfn = make_gsl_function(&mut lp,
                                              self.lower_bound,
                                              self.lower_bound + 1.0)?;
            bindings::gsl_integration_qagiu(&mut gslfn.function,
                                            self.lower_bound,
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

impl Integrator for QAGIL {
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
            let mut gslfn = make_gsl_function(&mut lp,
                                              self.upper_bound - 1.0,
                                              self.upper_bound)?;
            bindings::gsl_integration_qagil(&mut gslfn.function,
                                            self.upper_bound,
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
