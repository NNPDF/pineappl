#![allow(missing_docs)]

use pkg_config::Config;
use std::process::Command;

fn main() {
    let version = String::from_utf8(
        Command::new("fnlo-tk-config")
            .arg("--version")
            .output()
            .expect("did not find `fnlo-tk-config`, please install fastNLO")
            .stdout,
    )
    .unwrap();

    let tested_versions = ["2.5.0_2826"];

    if !tested_versions
        .iter()
        .any(|&tested| tested == version.trim())
    {
        println!(
            "cargo:warning=found fastNLO version {}, which has not been tested",
            version.trim()
        );
    }

    let fnlo_lib_path = String::from_utf8(
        Command::new("fnlo-tk-config")
            .arg("--libdir")
            .output()
            .expect("did not find `fnlo-tk-config`, please install fastNLO")
            .stdout,
    )
    .unwrap();

    println!("cargo:rustc-link-search={}", fnlo_lib_path.trim());

    let fnlo_include_path = String::from_utf8(
        Command::new("fnlo-tk-config")
            .arg("--incdir")
            .output()
            .expect("did not find `fnlo-tk-config`, please install fastNLO")
            .stdout,
    )
    .unwrap();

    let link_modifier = if cfg!(feature = "static") {
        "static="
    } else {
        ""
    };

    println!("cargo:rustc-link-lib={link_modifier}fastnlotoolkit");

    let lhapdf_include_paths = Config::new()
        .atleast_version("6")
        .statik(cfg!(feature = "static"))
        .cargo_metadata(false)
        .probe("lhapdf")
        .unwrap()
        .include_paths;

    cxx_build::bridge("src/lib.rs")
        .file("src/fastnlo.cpp")
        .include(fnlo_include_path.trim())
        .includes(lhapdf_include_paths)
        .std("c++17") // apparently not supported by MSVC, but fastNLO probably can't be compiled on Windows
        .compile("fnlo-bridge");

    println!("cargo:rerun-if-changed=src/lib.rs");
    println!("cargo:rerun-if-changed=src/fastnlo.cpp");
    println!("cargo:rerun-if-changed=src/fastnlo.hpp");
}
