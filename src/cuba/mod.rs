use std::{slice, vec};
use std::convert::From;
use std::os::raw::{c_int, c_longlong, c_void};

use super::traits::{IntegrandInput, IntegrandOutput};
use super::{IntegrationResult, Real};

mod cuhre;
pub use self::cuhre::Cuhre;

mod vegas;
pub use self::vegas::Vegas;

unsafe extern "C"
fn cuba_integrand<A, B, F>(ndim: *const c_int,
                           x: *const Real,
                           ncomp: *const c_int,
                           f: *mut Real,
                           userdata: *mut c_void) -> c_int
    where A: IntegrandInput,
          B: IntegrandOutput,
          F: FnMut(A) -> B
{
    let fnptr = userdata as *mut F;
    let fun: &mut F = &mut *fnptr;

    if (*ndim as usize) != A::input_size() {
        panic!("Internal error - wrong number of dimensions passed in");
    }

    let args = slice::from_raw_parts(x, *ndim as usize);
    let output = slice::from_raw_parts_mut(f, *ncomp as usize);

    fun(A::from_args(args)).into_args(output);
    0
}

#[derive(Copy, Clone, Debug)]
pub struct CubaIntegrationResult {
    pub value: Real,
    pub error: Real,
    pub prob: Real
}

#[derive(Clone, Debug)]
pub struct CubaIntegrationResults {
    pub nregions: c_int,
    pub neval: c_longlong,
    pub fail: c_int,
    pub results: Vec<CubaIntegrationResult>
}

pub struct CubaResultsIter {
    iter: vec::IntoIter<CubaIntegrationResult>
}

impl Iterator for CubaResultsIter {
    type Item = IntegrationResult;
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|CubaIntegrationResult { value, error, .. }| {
            IntegrationResult {
                value, error
            }
        })
    }
}

impl From<Vec<CubaIntegrationResult>> for CubaResultsIter {
    fn from(v: Vec<CubaIntegrationResult>) -> Self {
        CubaResultsIter {
            iter: v.into_iter()
        }
    }
}

impl super::traits::IntegrationResults for CubaIntegrationResults {
    type Iterator = CubaResultsIter;
    fn results(self) -> CubaResultsIter {
        From::from(self.results)
    }
}
