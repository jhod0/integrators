pub mod traits;

#[cfg(any(feature = "cuba", feature = "gsl"))]
mod bindings;

#[cfg(feature = "cuba")]
pub mod cuba;

#[cfg(feature = "cuba")]
pub use cuba::{Cuhre, Vegas};

#[cfg(feature = "gsl")]
pub mod gsl;

#[cfg(feature = "gsl")]
pub use gsl::{GSLIntegrationError, GSLResult,
              QNG, QAG, QAGRule, QAGS, QAGP};

#[cfg(test)]
mod test;

pub type Real = f64;

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
