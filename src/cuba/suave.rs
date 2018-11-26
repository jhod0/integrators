use std::{mem, ptr};
use std::os::raw::{c_int, c_longlong};
use ::bindings;
use ::traits::{IntegrandInput, IntegrandOutput};
use ::{Integrator, Real};

use super::{cuba_integrand, CubaError, CubaIntegrationResult, CubaIntegrationResults};

#[derive(Copy, Clone, Debug)]
pub struct Suave {
    mineval: usize,
    maxeval: usize,
    seed: usize,
    nnew: usize,
    nmin: usize,
    flatness: Real,
    flags: c_int,
}

impl Default for Suave {
    fn default() -> Self {
        Suave {
            mineval: 1,
            maxeval: c_longlong::max_value() as usize,
            seed: 0,
            nnew: 1000,
            nmin: 5,
            flatness: 25 as Real,
            flags: 0,
        }
    }
}

impl Suave {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_mineval(self, maxeval: usize) -> Self {
        Suave {
            maxeval, ..self
        }
    }

    pub fn with_maxeval(self, maxeval: usize) -> Self {
        Suave {
            maxeval, ..self
        }
    }

    pub fn with_seed(self, seed: usize) -> Self {
        Suave {
            seed, ..self
        }
    }

    pub fn with_nnew(self, nnew: usize) -> Self {
        Suave {
            nnew, ..self
        }
    }

    pub fn with_nmin(self, nmin: usize) -> Self {
        Suave {
            nmin, ..self
        }
    }

    pub fn with_flatness(self, flatness: Real) -> Self {
        Suave {
            flatness, ..self
        }
    }
}

impl Integrator for Suave {
    type Success = CubaIntegrationResults;
    type Failure = super::CubaError;
    fn integrate<A, B, F: FnMut(A) -> B>(&mut self, mut fun: F, epsrel: Real, epsabs: Real) -> Result<Self::Success, Self::Failure>
        where A: IntegrandInput,
              B: IntegrandOutput
    {
        // Using cuba's parallelization via fork() would deeply break Rust's
        // concurrency model and safety guarantees. So, we'll turn it off.
        unsafe { bindings::cubacores(0, 0) };

        let (ndim, ncomp) = {
            let inputs = A::input_size();
            let outputs = fun(A::from_args(&vec![0.5; inputs][..])).output_size();
            (inputs, outputs)
        };

        let mut nregions = 0;
        let mut neval = 0;
        let mut fail = 0;
        let (mut value, mut error, mut prob) =
                (vec![0.0; ncomp], vec![0.0; ncomp], vec![0.0; ncomp]);

        unsafe {
            bindings::llSuave(ndim as c_int, ncomp as c_int,
                              Some(cuba_integrand::<A, B, F>), mem::transmute(&mut fun),
                              1 /* nvec */,
                              epsrel,
                              epsabs,
                              self.flags,
                              self.seed as c_int,
                              self.mineval as c_longlong,
                              self.maxeval as c_longlong,
                              self.nnew as c_longlong,
                              self.nmin as c_longlong,
                              self.flatness,
                              // statefile
                              ptr::null(),
                              // spin
                              ptr::null_mut(),
                              &mut nregions,
                              &mut neval,
                              &mut fail,
                              value.as_mut_ptr(),
                              error.as_mut_ptr(),
                              prob.as_mut_ptr());
        }

        if fail == 0 {
            Ok(CubaIntegrationResults {
                nregions: Some(nregions), neval,
                results: value.iter().zip(error.iter()).zip(prob.iter())
                              .map(|((&value, &error), &prob)|
                                     CubaIntegrationResult {
                                         value, error, prob
                                     })
                              .collect()
            })
        } else if fail == -1 {
            // `baddim`
            Err(CubaError::BadDim("suave", ndim))
        } else if fail == -2 {
            // `badcomp`
            Err(CubaError::BadComp("suave", ncomp))
        } else if fail == 1 {
            Err(CubaError::DidNotConverge(CubaIntegrationResults {
                nregions: None, neval,
                results: value.iter().zip(error.iter()).zip(prob.iter())
                              .map(|((&value, &error), &prob)|
                                     CubaIntegrationResult {
                                         value, error, prob
                                     })
                              .collect()
            }))
        } else {
            unreachable!("Suave returned invalid failure code: {}", fail)
        }
    }
}
