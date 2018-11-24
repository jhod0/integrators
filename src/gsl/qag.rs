use std::convert::Into;
use std::os::raw::c_int;

use ::bindings;
use ::{Integrator, Real};
use ::traits::{IntegrandInput, IntegrandOutput};

use super::{make_gsl_function, GSLIntegrationResult, GSLIntegrationWorkspace};

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

pub struct QAG {
    range_low: Real,
    range_high: Real,
    wkspc: GSLIntegrationWorkspace,
    rule: QAGRule
}

impl QAG {
    pub fn new(nintervals: usize) -> Self {
        QAG {
            range_low: 0.0,
            range_high: 1.0,
            rule: QAGRule::Gauss61,
            wkspc: GSLIntegrationWorkspace::new(nintervals)
        }
    }

    pub fn with_nintervals(self, nintervals: usize) -> Self {
        QAG {
            wkspc: GSLIntegrationWorkspace::new(nintervals),
            ..self
        }
    }

    pub fn with_range(self, range_low: Real, range_high: Real) -> Self {
        QAG { range_low, range_high, ..self }
    }

    pub fn with_rule(self, rule: QAGRule) -> Self {
        QAG { rule, ..self }
    }
}

impl Integrator for QAG {
    type Success = GSLIntegrationResult;
    type Failure = ();
    fn integrate<A, B, F: FnMut(A) -> B>(&mut self, mut fun: F, epsrel: Real, epsabs: Real) -> Result<Self::Success, Self::Failure>
        where A: IntegrandInput,
              B: IntegrandOutput
    {
        let mut gslfn = make_gsl_function(&mut fun);
        let mut value: Real = 0.0;
        let mut error: Real = 0.0;

        let retcode = unsafe {
            bindings::gsl_integration_qag(&mut gslfn.function,
                                          self.range_low, self.range_high,
                                          epsabs, epsrel,
                                          self.wkspc.nintervals,
                                          self.rule.into(),
                                          self.wkspc.wkspc,
                                          &mut value,
                                          &mut error)
        };

        if retcode != bindings::GSL_SUCCESS {
            panic!("GSL failed");
        } else {
            Ok(GSLIntegrationResult {
                value, error
            })
        }
    }
}
