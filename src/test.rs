use super::{Integrator, Real};
#[cfg(feature = "cuba")]
use super::cuba::{Cuhre, CubaError, Vegas};

#[test]
#[cfg(feature = "cuba")]
fn test_simple_integration() {
    let mut cuhre = Cuhre::new(1000000);
    let mut vegas = Vegas::default().with_maxeval(1000000);

    let a = cuhre.integrate(|a: Real| (a * a),
                            1e-4, 1e-12);
    assert_eq!(a, Err(CubaError::BadDim("cuhre", 1)));

    let b = vegas.integrate(|a: Real| (a * a),
                            1e-4, 1e-12);
    assert!(b.is_ok());
}
