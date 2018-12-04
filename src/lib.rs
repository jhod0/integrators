pub mod traits;
pub mod ffi;

#[cfg(any(feature = "cuba", feature = "gsl"))]
mod bindings;

#[cfg(feature = "cuba")]
pub mod cuba;

#[cfg(feature = "gsl")]
pub mod gsl;

#[cfg(test)]
mod test;

pub type Real = f64;
pub type Real2 = (Real, Real);
pub type Real3 = (Real, Real, Real);
pub type Real4 = (Real, Real, Real, Real);
pub type Real5 = (Real, Real, Real, Real, Real);
pub type Real6 = (Real, Real, Real, Real, Real, Real);
pub type Real7 = (Real, Real, Real, Real, Real, Real, Real);
pub type Real8 = (Real, Real, Real, Real, Real, Real, Real, Real);

pub use traits::{Integrator, IntegrandInput, IntegrandOutput};

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct IntegrationResult {
    pub value: Real,
    pub error: Real,
}

pub struct IntegrationResultIter {
    val: Option<IntegrationResult>
}

impl Iterator for IntegrationResultIter {
    type Item = IntegrationResult;
    fn next(&mut self) -> Option<Self::Item> {
        let mut output: Option<Self::Item> = None;
        self.val = match self.val {
            Some(item) => {
                output = Some(item);
                None
            },
            None => None
        };
        output
    }
}

impl self::traits::IntegrationResults for IntegrationResult {
    type Iterator = IntegrationResultIter;
    fn results(self) -> Self::Iterator {
        IntegrationResultIter {
            val: Some(self)
        }
    }
}
