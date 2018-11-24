use super::{Integrator, Real};
#[cfg(feature = "cuba")]
use super::{Cuhre, Vegas};

#[test]
#[cfg(feature = "cuba")]
fn test_simple_integration() {
    let mut cuhre = Cuhre::new(1000000);
    let mut vegas = Vegas::default().with_maxeval(1000000);

    let a = cuhre.integrate(|a: Real| (a * a),
                            1e-4, 1e-12);
    println!("{:?}", a);

    let b = vegas.integrate(|a: Real| (a * a),
                            1e-4, 1e-12);
    println!("{:?}", b);

    let c = vegas.integrate(|a: Real| {
        if (a >= 0.4) && (a <= 0.6) {
            0.0 as Real
        } else {
            panic!("Boom")
        }},
        1e-4, 1e-12);
    println!("{:?}", c);
}
