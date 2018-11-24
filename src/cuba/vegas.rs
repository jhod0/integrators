use std::{mem, ptr};
use std::os::raw::{c_int, c_longlong};
use ::bindings;
use ::traits::{IntegrandInput, IntegrandOutput};
use ::{Integrator, Real};

use super::{cuba_integrand, CubaIntegrationResult, CubaIntegrationResults};

#[derive(Copy, Clone, Debug)]
pub struct Vegas {
    mineval: usize,
    maxeval: usize,
    seed: usize,
    nstart: usize,
    nincrease: usize,
    nbatch: usize,
    flags: c_int,
}

impl Default for Vegas {
    fn default() -> Self {
        Vegas {
            mineval: 1,
            maxeval: c_longlong::max_value() as usize,
            seed: 0,
            nstart: 1000,
            nincrease: 500,
            nbatch: 1000,
            flags: 0
        }
    }
}

impl Vegas {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_mineval(self, maxeval: usize) -> Self {
        Vegas {
            maxeval, ..self
        }
    }

    pub fn with_maxeval(self, maxeval: usize) -> Self {
        Vegas {
            maxeval, ..self
        }
    }

    pub fn with_seed(self, seed: usize) -> Self {
        Vegas {
            seed, ..self
        }
    }

    pub fn with_nstart(self, nstart: usize) -> Self {
        Vegas {
            nstart, ..self
        }
    }

    pub fn with_nincrease(self, nincrease: usize) -> Self {
        Vegas {
            nincrease, ..self
        }
    }

    pub fn with_nbatch(self, nbatch: usize) -> Self {
        Vegas {
            nbatch, ..self
        }
    }
}

impl Integrator for Vegas {
    type Success = CubaIntegrationResults;
    type Failure = super::CubaError;
    fn integrate<A, B, F: FnMut(A) -> B>(&mut self, mut fun: F, epsrel: Real, epsabs: Real) -> Result<Self::Success, Self::Failure>
        where A: IntegrandInput,
              B: IntegrandOutput
    {
        let (ndim, ncomp) = {
            let inputs = A::input_size();
            let outputs = fun(A::from_args(&vec![0.5; inputs][..])).output_size();
            (inputs, outputs)
        };

        let mut neval = 0;
        let mut fail = 0;
        let (mut value, mut error, mut prob) =
                (vec![0.0; ncomp], vec![0.0; ncomp], vec![0.0; ncomp]);

        unsafe {
            bindings::llVegas(ndim as c_int, ncomp as c_int,
                              Some(cuba_integrand::<A, B, F>), mem::transmute(&mut fun),
                              1 /* nvec */,
                              epsrel,
                              epsabs,
                              self.flags,
                              self.seed as c_int,
                              self.mineval as c_longlong,
                              self.maxeval as c_longlong,
                              self.nstart as c_longlong,
                              self.nincrease as c_longlong,
                              self.nbatch as c_longlong,
                              0 /* gridno */,
                              // statefile
                              ptr::null(),
                              // spin
                              ptr::null_mut(),
                              &mut neval,
                              &mut fail,
                              value.as_mut_ptr(),
                              error.as_mut_ptr(),
                              prob.as_mut_ptr());
        }

        Ok(CubaIntegrationResults {
            nregions: -1, neval, fail,
            results: value.iter().zip(error.iter()).zip(prob.iter())
                          .map(|((&value, &error), &prob)|
                                 CubaIntegrationResult {
                                     value, error, prob
                                 })
                          .collect()
        })
    }
}
