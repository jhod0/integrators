use std::{mem, ptr};
use std::os::raw::{c_int, c_longlong};

use ::bindings;
use ::ffi::LandingPad;
use ::traits::{IntegrandInput, IntegrandOutput};
use ::{Integrator, Real};

use super::{cuba_integrand, CubaError, CubaIntegrationResult, CubaIntegrationResults};

#[derive(Copy, Clone, Debug)]
pub struct Cuhre {
    pub mineval: usize,
    pub maxeval: usize,
    key: Option<u16>,
}

impl Cuhre {
    pub fn new(maxeval: usize) -> Self {
        Cuhre {
            mineval: 1, maxeval, key: None
        }
    }

    pub fn with_mineval(self, mineval: usize) -> Self {
        Cuhre {
            mineval, ..self
        }
    }

    pub fn with_maxeval(self, maxeval: usize) -> Self {
        Cuhre {
            maxeval, ..self
        }
    }

    pub fn with_key(self, key: Option<u16>) -> Option<Self> {
        if key.map(|k| [7, 9, 11, 13].contains(&k)).unwrap_or(true) {
            None
        } else {
            Some(Cuhre {
                key, ..self
            })
        }
    }
}

impl Integrator for Cuhre {
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

        let key = match (self.key, ndim) {
            (Some(key), _) => key,
            (_, 1) | (_, 2) => 13,
            (_, 3) => 11,
            _ => 9
        };

        assert!([7, 9, 11, 13].contains(&key));

        let mut lp = LandingPad::new(fun);
        unsafe {
            bindings::llCuhre(ndim as c_int, ncomp as c_int,
                              Some(cuba_integrand::<A, B, F>), mem::transmute(&mut lp),
                              1 /* nvec */,
                              epsrel,
                              epsabs,
                              0 /* flags */,
                              self.mineval as c_longlong,
                              self.maxeval as c_longlong,
                              key as c_int,
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
        lp.maybe_resume_unwind();

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
            Err(CubaError::BadDim("cuhre", ndim))
        } else if fail == -2 {
            // `badcomp`
            Err(CubaError::BadComp("cuhre", ncomp))
        } else if fail == 1 {
            Err(CubaError::DidNotConverge(CubaIntegrationResults {
                nregions: Some(nregions), neval,
                results: value.iter().zip(error.iter()).zip(prob.iter())
                              .map(|((&value, &error), &prob)|
                                     CubaIntegrationResult {
                                         value, error, prob
                                     })
                              .collect()
            }))
        } else {
            unreachable!("panic should have been propogated into Rust caller")
        }
    }
}
