use std::convert::Into;
use std::os::raw::c_int;

use ::bindings;
use ::{IntegrationResult, Integrator, Real};
use ::ffi::LandingPad;
use ::traits::{IntegrandInput, IntegrandOutput};

use super::{make_gsl_function, GSLIntegrationError, GSLIntegrationWorkspace};

/// Quadrature rule to apply for QAG integration. Rules are supported for 15,
/// 21, 31, 41, 51, 61 points.
#[derive(Debug, Hash, Copy, Clone, PartialEq, Eq)]
pub enum QAGRule {
    Gauss15,
    Gauss21,
    Gauss31,
    Gauss41,
    Gauss51,
    Gauss61
}

impl Into<c_int> for QAGRule {
    fn into(self) -> c_int {
        use self::QAGRule::*;
        match self {
            Gauss15 => bindings::GSL_INTEG_GAUSS15 as c_int,
            Gauss21 => bindings::GSL_INTEG_GAUSS21 as c_int,
            Gauss31 => bindings::GSL_INTEG_GAUSS31 as c_int,
            Gauss41 => bindings::GSL_INTEG_GAUSS41 as c_int,
            Gauss51 => bindings::GSL_INTEG_GAUSS51 as c_int,
            Gauss61 => bindings::GSL_INTEG_GAUSS61 as c_int,
        }
    }
}

/// Quadrature Adaptive General integration. Iteratively applies quadrature
/// rules to sub-regions of a function until it converges to the requested
/// uncertainty. This method will work reasonably well for most well-behaved
/// functions.
///
/// ```
/// use integrators::{Integrator, Real};
/// let mut qag = integrators::QAG::new(1000);
///
/// let res1 = qag.integrate(|a: Real| a * a, 1e-6, 1e-10)
///               .unwrap();
///
/// assert!((res1.value - 3f64.recip()).abs() < res1.error);
///
/// let res2 = qag.with_range(3.0, 10.0)
///               .integrate(|a: Real| a * a, 1e-6, 1e-10)
///               .unwrap();
/// assert!((res2.value - (1000f64 - 27f64) / 3.).abs() < res2.error);
/// ```
#[derive(Debug, Clone)]
pub struct QAG {
    range_low: Real,
    range_high: Real,
    wkspc: GSLIntegrationWorkspace,
    rule: QAGRule
}

impl QAG {
    /// Creates a new QAG with enough memory for `nintervals` subintervals.
    /// This will create a QAG to integrate the range [0, 1], and using the
    /// 61-point (highest-order) Gauss-Kronrod rule. To change the integration
    /// bounds, see `with_range`, and to change the quadrature rule, see
    /// `with_rule`.
    pub fn new(nintervals: usize) -> Self {
        QAG {
            range_low: 0.0,
            range_high: 1.0,
            rule: QAGRule::Gauss61,
            wkspc: GSLIntegrationWorkspace::new(nintervals)
        }
    }

    /// Discards the old workspace and allocates a new one with enough memory
    /// for `nintervals` subintervals.
    pub fn with_nintervals(self, nintervals: usize) -> Self {
        QAG {
            wkspc: GSLIntegrationWorkspace::new(nintervals),
            ..self
        }
    }

    /// Use a different integration range. (Default = [0, 1])
    pub fn with_range(self, range_low: Real, range_high: Real) -> Self {
        QAG { range_low, range_high, ..self }
    }


    /// Use a specific quadrature rule. (Default = `QAGRule::Gauss61`)
    pub fn with_rule(self, rule: QAGRule) -> Self {
        QAG { rule, ..self }
    }
}

impl Integrator for QAG {
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
            bindings::gsl_integration_qag(&mut gslfn.function,
                                          self.range_low, self.range_high,
                                          epsabs, epsrel,
                                          self.wkspc.nintervals,
                                          self.rule.into(),
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
