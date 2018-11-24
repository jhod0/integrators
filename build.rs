extern crate bindgen;

use std::env;
use std::path::PathBuf;

fn main() {
    println!("cargo:rustc-link-lib=cuba");

    let mut bindings = bindgen::Builder::default();
    let mut any_bindings = false;

    if env::var("CARGO_FEATURE_GSL").is_ok() {
        bindings = bindings.header("gsl_wrapper.h");
        any_bindings = true;
    }
    if env::var("CARGO_FEATURE_CUBA").is_ok() {
        bindings = bindings.header("cuba_wrapper.h");
        any_bindings = true;
    }

    if any_bindings {
        let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
        bindings.blacklist_item("FP_NORMAL")
                .blacklist_item("FP_SUBNORMAL")
                .blacklist_item("FP_ZERO")
                .blacklist_item("FP_INFINITE")
                .blacklist_item("FP_NAN")
                .generate().unwrap()
                .write_to_file(out_path.join("integrand_bindings.rs"))
                .expect("Error writing bindings");
    }
}
