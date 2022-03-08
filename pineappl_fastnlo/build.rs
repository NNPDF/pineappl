use std::process::Command;

fn main() {
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

    println!("cargo:rustc-link-lib=fastnlotoolkit");

    cxx_build::bridge("src/lib.rs")
        .file("src/fastnlo.cpp")
        .include(fnlo_include_path.trim())
        .compile("fnlo-bridge");

    println!("cargo:rerun-if-changed=src/lib.rs");
    println!("cargo:rerun-if-changed=src/fastnlo.cpp");
    println!("cargo:rerun-if-changed=src/fastnlo.hpp");
}
