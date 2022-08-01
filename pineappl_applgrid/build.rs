use std::process::Command;

fn main() {
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

    cxx_build::bridge("src/lib.rs")
        .file("src/applgrid.cpp")
        .include(include_path.trim())
        .compile("appl-bridge");

    println!("cargo:rerun-if-changed=src/lib.rs");
    println!("cargo:rerun-if-changed=src/applgrid.cpp");
    println!("cargo:rerun-if-changed=src/applgrid.hpp");
    println!("cargo:rerun-if-changed=src/calculation.hpp");
}
