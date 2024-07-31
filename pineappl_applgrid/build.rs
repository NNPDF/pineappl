#![allow(missing_docs)]

use cc::Build;
use pkg_config::Config;
use std::env;
use std::path::Path;
use std::process::Command;

fn conditional_std<'a>(build: &'a mut Build, std: Option<&str>) -> &'a mut Build {
    if let Some(std) = std {
        build.std(std)
    } else {
        build
    }
}

fn main() {
    let version = String::from_utf8(
        Command::new("applgrid-config")
            .arg("--version")
            .output()
            .expect("did not find `applgrid-config`, please install APPLgrid")
            .stdout,
    )
    .unwrap();

    let tested_versions = [
        "1.6.27", "1.6.28", "1.6.29", "1.6.30", "1.6.31", "1.6.32", "1.6.35", "1.6.36",
    ];

    if !tested_versions
        .iter()
        .any(|&tested| tested == version.trim())
    {
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

    let appl_igrid_dir = env::var("APPL_IGRID_DIR").unwrap_or_else(|_| {
        Path::new(
            &String::from_utf8(
                Command::new("applgrid-config")
                    .arg("--incdir")
                    .output()
                    .expect("did not find `applgrid-config`, please install APPLgrid")
                    .stdout,
            )
            .unwrap(),
        )
        .join("appl_grid")
        .to_str()
        .unwrap()
        .to_owned()
    });

    let cxx_flags: Vec<_> = String::from_utf8(
        Command::new("applgrid-config")
            .arg("--cxxflags")
            .output()
            .expect("did not find `applgrid-config`, please install APPLgrid")
            .stdout,
    )
    .unwrap()
    .split_ascii_whitespace()
    .map(ToOwned::to_owned)
    .collect();

    let include_dirs: Vec<_> = cxx_flags
        .iter()
        .filter_map(|token| token.strip_prefix("-I"))
        .collect();

    let std = cxx_flags
        .iter()
        .filter_map(|token| token.strip_prefix("-std="))
        .last();

    let libs = String::from_utf8(
        Command::new("applgrid-config")
            .arg("--ldflags")
            .output()
            .expect("did not find `applgrid-config`, please install APPLgrid")
            .stdout,
    )
    .unwrap();

    let link_modifier = if cfg!(feature = "static") {
        // for some reason `libz.a` isn't found, although `libz.so` is
        let zlib_link_paths = Config::new()
            .cargo_metadata(false)
            .statik(true)
            .probe("zlib")
            .unwrap()
            .link_paths;

        for link_path in zlib_link_paths {
            println!("cargo:rustc-link-search={}", link_path.to_str().unwrap());
        }

        "static="
    } else {
        ""
    };

    for lib in libs
        .split_whitespace()
        .filter_map(|token| token.strip_prefix("-l"))
    {
        match lib {
            // we can't link gfortran statically - to avoid it compile APPLgrid without HOPPET
            "gfortran" => println!("cargo:rustc-link-lib={lib}"),
            _ => println!("cargo:rustc-link-lib={link_modifier}{lib}"),
        }
    }

    conditional_std(
        Build::new()
            .cpp(true)
            .file("src/check_appl_igrid.cpp")
            .includes(&include_dirs)
            .include(&appl_igrid_dir),
        std,
    )
    .try_compile("appl_igrid")
    .expect(
        "could not find file `appl_igrid.h`, please set the environment variable \
                `APPL_IGRID_DIR` to the directory containing it",
    );

    println!("cargo:rerun-if-env-changed=APPL_IGRID_DIR");

    conditional_std(
        cxx_build::bridge("src/lib.rs")
            .file("src/applgrid.cpp")
            .includes(&include_dirs)
            .include(appl_igrid_dir),
        std,
    )
    .compile("appl-bridge");

    println!("cargo:rerun-if-changed=src/lib.rs");
    println!("cargo:rerun-if-changed=src/applgrid.cpp");
    println!("cargo:rerun-if-changed=src/applgrid.hpp");
    println!("cargo:rerun-if-changed=src/helpers.hpp");
    println!("cargo:rerun-if-changed=src/check_appl_igrid.cpp");
}
