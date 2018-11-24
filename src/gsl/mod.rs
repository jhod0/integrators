use std::{marker, mem};
use std::os::raw::c_void;

use super::bindings;
use super::traits::{IntegrandInput, IntegrandOutput};
use super::Real;

mod qag;
pub use self::qag::{QAG, QAGRule};

mod qng;
pub use self::qng::QNG;

unsafe extern "C"
fn gsl_integrand_fn<A, B, F>(x: Real, params: *mut c_void) -> Real
    where A: IntegrandInput,
          B: IntegrandOutput,
          F: FnMut(A) -> B
{
    let fnptr = params as *mut F;
    let fun: &mut F = &mut *fnptr;

    if A::input_size() != 1 {
        panic!("integrand given to GSL integrator demands >1 input");
    }

    let mut output: [Real; 1] = [0.0];
    fun(A::from_args(&[x])).into_args(&mut output);
    output[0]
}

struct GSLFunction<'a> {
    function: bindings::gsl_function,
    lifetime: marker::PhantomData<&'a ()>
}

fn make_gsl_function<'a, A, B, F>(fun: &'a mut F) -> GSLFunction<'a>
    where A: IntegrandInput,
          B: IntegrandOutput,
          F: FnMut(A) -> B
{
    GSLFunction {
        function: bindings::gsl_function {
            function: Some(gsl_integrand_fn::<A, B, F>),
            params: unsafe { mem::transmute(fun) }
        },
        lifetime: marker::PhantomData
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct GSLIntegrationResult {
    pub value: Real,
    pub error: Real
}

struct GSLIntegrationWorkspace {
    pub(crate) nintervals: usize,
    wkspc: *mut bindings::gsl_integration_workspace
}

impl GSLIntegrationWorkspace {
    pub(crate) fn new(n: usize) -> Self {
        GSLIntegrationWorkspace {
            wkspc: unsafe {
                bindings::gsl_integration_workspace_alloc(n)
            },
            nintervals: n
        }
    }
}

impl Drop for GSLIntegrationWorkspace {
    fn drop(&mut self) {
        unsafe {
            bindings::gsl_integration_workspace_free(self.wkspc)
        }
    }
}
