use std::error;
use super::{Real, IntegrationResult};

pub trait Integrator {
    type Success: IntegrationResults;
    type Failure: error::Error;
    fn integrate<A, B, F: FnMut(A) -> B>(&mut self, fun: F, epsrel: Real, epsabs: Real) -> Result<Self::Success, Self::Failure>
        where A: IntegrandInput,
              B: IntegrandOutput;
}

pub trait IntegrandInput {
    fn input_size() -> usize;
    fn from_args(&[Real]) -> Self;
}

pub trait IntegrandOutput {
    fn output_size(&self) -> usize;
    fn into_args(&self, &mut [Real]);
}

pub trait IntegrationResults {
    type Iterator: Iterator<Item=IntegrationResult>;
    fn results(self) -> Self::Iterator;
}

impl IntegrandOutput for Vec<Real> {
    fn output_size(&self) -> usize {
        self.len()
    }

    fn into_args(&self, output: &mut [Real]) {
        if self.len() != output.len() {
            panic!("Integrand returned vector of wrong length: expected {}, got {}",
                   output.len(), self.len());
        }
        for i in 0..self.len() {
            output[i] = self[i];
        }
    }
}

macro_rules! impl_integrand_traits {
    ($ty:ty, $size:expr, $convert:expr, $result:expr) => {
        impl IntegrandInput for $ty {
            fn input_size() -> usize {
                $size
            }

            fn from_args(args: &[Real]) -> Self {
                assert!(args.len() == Self::input_size());
                $convert(args)
            }
        }

        impl IntegrandOutput for $ty {
            fn output_size(&self) -> usize {
                $size
            }

            fn into_args(&self, args: &mut [Real]) {
                assert!(args.len() == self.output_size());
                $result(self, args)
            }
        }
    }
}

impl_integrand_traits!(Real, 1,
                       |args: &[Real]| { args[0] },
                       |this: &Real, args: &mut [Real]| {
                           args[0] = *this;
                       });
impl_integrand_traits!((Real, Real), 2,
                       |args: &[Real]| { (args[0], args[1]) },
                       |this: &(Real, Real), args: &mut [Real]| {
                           args[0] = this.0;
                           args[1] = this.1;
                       });
impl_integrand_traits!((Real, Real, Real), 3,
                       |args: &[Real]| { (args[0], args[1], args[2]) },
                       |this: &(Real, Real, Real), args: &mut [Real]| {
                           args[0] = this.0;
                           args[1] = this.1;
                           args[2] = this.2;
                       });
impl_integrand_traits!((Real, Real, Real, Real), 4,
                       |args: &[Real]| { (args[0], args[1], args[2], args[3]) },
                       |this: &(Real, Real, Real, Real), args: &mut [Real]| {
                           args[0] = this.0;
                           args[1] = this.1;
                           args[2] = this.2;
                           args[3] = this.3;
                       });
impl_integrand_traits!((Real, Real, Real, Real, Real), 5,
                       |args: &[Real]| { (args[0], args[1], args[2], args[3], args[4]) },
                       |this: &(Real, Real, Real, Real, Real), args: &mut [Real]| {
                           args[0] = this.0;
                           args[1] = this.1;
                           args[2] = this.2;
                           args[3] = this.3;
                           args[4] = this.4;
                       });
impl_integrand_traits!((Real, Real, Real, Real, Real, Real), 6,
                       |args: &[Real]| { (args[0], args[1], args[2], args[3], args[4], args[5]) },
                       |this: &(Real, Real, Real, Real, Real, Real), args: &mut [Real]| {
                           args[0] = this.0;
                           args[1] = this.1;
                           args[2] = this.2;
                           args[3] = this.3;
                           args[4] = this.4;
                           args[5] = this.5;
                       });
impl_integrand_traits!((Real, Real, Real, Real, Real, Real, Real), 7,
                       |args: &[Real]| { (args[0], args[1], args[2], args[3], args[4], args[5], args[6]) },
                       |this: &(Real, Real, Real, Real, Real, Real, Real), args: &mut [Real]| {
                           args[0] = this.0;
                           args[1] = this.1;
                           args[2] = this.2;
                           args[3] = this.3;
                           args[4] = this.4;
                           args[5] = this.5;
                           args[6] = this.6;
                       });
impl_integrand_traits!((Real, Real, Real, Real, Real, Real, Real, Real), 7,
                       |args: &[Real]| { (args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7]) },
                       |this: &(Real, Real, Real, Real, Real, Real, Real, Real), args: &mut [Real]| {
                           args[0] = this.0;
                           args[1] = this.1;
                           args[2] = this.2;
                           args[3] = this.3;
                           args[4] = this.4;
                           args[5] = this.5;
                           args[6] = this.6;
                           args[7] = this.7;
                       });


#[cfg(test)]
mod test_traits {
    use super::{Real, IntegrandInput, IntegrandOutput};

    #[test]
    fn test_from_into_traits() {
        let a: Real = Real::from_args(&[2.0]);
        assert_eq!(a, 2.0);
        let mut b: [Real; 1] = [0.0];
        a.into_args(&mut b);
        assert_eq!(b[0], a);

        let c: (Real, Real) = <(Real, Real)>::from_args(&[1.0, 2.0]);
        assert_eq!(c, (1.0, 2.0));
        let mut d: [Real; 2] = [0.0, 0.0];
        c.into_args(&mut d);
        assert_eq!((d[0], d[1]), c);
    }

    #[test]
    #[should_panic]
    fn test_from_failure() {
        let _a: Real = Real::from_args(&[2.0, 3.0]);
    }

    #[test]
    #[should_panic]
    fn test_into_failure() {
        let args: &'static [Real] = &[0.1, 0.2, 0.3, 0.4];
        let _c: (Real, Real, Real) = <(Real, Real, Real)>::from_args(args);
    }

    #[test]
    fn test_vec_traits() {
        let v: Vec<Real> = vec![0.5; 10];
        let mut args: [Real; 10] = [0.0; 10];
        v.into_args(&mut args);
        for (&a, &b) in v.iter().zip(args.iter()) {
            assert_eq!(a, b)
        }
    }

    #[test]
    #[should_panic(expected = "Integrand returned vector of wrong length: expected 10, got 11")]
    fn test_vec_traits_failure() {
        let v: Vec<Real> = vec![0.5; 11];
        let mut args: [Real; 10] = [0.0; 10];
        v.into_args(&mut args);
    }
}
