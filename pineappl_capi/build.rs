//! A build script to install the OOP C++ interface to `PineAPPL`.

use std::env;
use std::fs;
use std::path::PathBuf;

fn main() {
    println!("cargo:rerun-if-changed=src/include/PineAPPL.hpp");

    if let Ok(prefix) = env::var("PINEAPPL_C_INSTALL_PREFIX") {
        let prefix_path = PathBuf::from(prefix);
        let include_path = prefix_path.join("include").join("pineappl_capi");

        fs::create_dir_all(&include_path).expect("Failed to create include directory.");

        let source_header = PathBuf::from("src/include/PineAPPL.hpp");
        let dest_header = include_path.join("PineAPPL.hpp");

        fs::copy(&source_header, &dest_header).expect("Failed to copy header file.");
    }
}
