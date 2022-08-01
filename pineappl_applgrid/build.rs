use cc::Build;
use std::env;
use std::path::Path;
use std::process::Command;

fn main() {
    let version = String::from_utf8(
        Command::new("applgrid-config")
            .arg("--version")
            .output()
            .expect("did not find `applgrid-config`, please install APPLgrid")
            .stdout,
    )
    .unwrap();

    if version.trim() != "1.6.27" {
        println!(
            "cargo:warning=found APPLgrid version {}, which has not been tested",
            version.trim()
        );
    }

    let lib_path = String::from_utf8(
        Command::new("applgrid-config")
            .arg("--libdir")
            .output()
            .expect("did not find `applgrid-config`, please install APPLgrid")
            .stdout,
    )
    .unwrap();

    println!("cargo:rustc-link-search={}", lib_path.trim());

    let include_path = String::from_utf8(
        Command::new("applgrid-config")
            .arg("--incdir")
            .output()
            .expect("did not find `applgrid-config`, please install APPLgrid")
            .stdout,
    )
    .unwrap();

    let appl_igrid_dir = if let Ok(dir) = env::var("APPL_IGRID_DIR") {
        dir
    } else {
        Path::new(&include_path)
            .join("appl_grid")
            .to_str()
            .unwrap()
            .to_owned()
    };

    let libs = String::from_utf8(
        Command::new("applgrid-config")
            .arg("--ldflags")
            .output()
            .expect("did not find `applgrid-config`, please install APPLgrid")
            .stdout,
    )
    .unwrap();

    for lib in libs
        .split_whitespace()
        .filter_map(|token| token.strip_prefix("-l"))
    {
        println!("cargo:rustc-link-lib={}", lib);
    }

    Build::new()
        .file("src/check_appl_igrid.cpp")
        .include(&appl_igrid_dir)
        .try_compile("appl_igrid")
        .expect(
            "could not find file `appl_igrid.h`, please set the environment variable \
                `APPL_IGRID_DIR` to the directory containing it",
        );

    println!("cargo:rerun-if-env-changed=APPL_IGRID_DIR");

    cxx_build::bridge("src/lib.rs")
        .file("src/applgrid.cpp")
        .include(include_path.trim())
        .include(appl_igrid_dir)
        .compile("appl-bridge");

    println!("cargo:rerun-if-changed=src/lib.rs");
    println!("cargo:rerun-if-changed=src/applgrid.cpp");
    println!("cargo:rerun-if-changed=src/applgrid.hpp");
    println!("cargo:rerun-if-changed=src/calculation.hpp");
    println!("cargo:rerun-if-changed=src/check_appl_igrid_dir.cpp");
}
