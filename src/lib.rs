pub mod traits;

#[cfg(any(feature = "cuba", feature = "gsl"))]
mod bindings;

#[cfg(feature = "cuba")]
mod cuba;

#[cfg(feature = "cuba")]
pub use cuba::{Cuhre, Vegas};

#[cfg(feature = "gsl")]
mod gsl;

#[cfg(feature = "gsl")]
pub use gsl::{GSLIntegrationResult, QAG, QAGRule, QNG};

#[cfg(test)]
mod test;

pub type Real = f64;

pub use traits::{Integrator, IntegrandInput, IntegrandOutput};

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct IntegrationResult {
    pub value: Real,
    pub error: Real,
}
