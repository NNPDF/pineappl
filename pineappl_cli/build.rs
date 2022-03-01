#[cfg(feature = "fastnlo")]
fn main() {
    use std::process::Command;

    let mut bridge = cxx_build::bridge("src/import.rs");

    // PineAPPL CAPI
    let pineappl_capi = pkg_config::Config::new()
        .atleast_version("0.5.0")
        .probe("pineappl_capi")
        .expect("PineAPPL's C API not found, please install it");

    for include_path in pineappl_capi.include_paths {
        bridge.include(include_path);
    }

    for lib_path in pineappl_capi.link_paths {
        println!("cargo:rustc-link-search={}", lib_path.to_str().unwrap());
    }

    for lib in pineappl_capi.libs {
        println!("cargo:rustc-link-lib={}", lib);
    }

    // fastNLO
    let fnlo_include_path = String::from_utf8(
        Command::new("fnlo-tk-config")
            .arg("--incdir")
            .output()
            .expect("fastNLO not found, please install it")
            .stdout,
    )
    .unwrap();

    bridge.include(fnlo_include_path.trim());

    let fnlo_lib_path = String::from_utf8(
        Command::new("fnlo-tk-config")
            .arg("--libdir")
            .output()
            .unwrap()
            .stdout,
    )
    .unwrap();

    println!("cargo:rustc-link-search={}", fnlo_lib_path.trim());

    // TODO: why do I have to link statically?
    println!("cargo:rustc-link-lib=fastnlotoolkit");

    // bridging code
    bridge.file("src/import/fastnlo.cpp").compile("cxx-bridge");

    println!("cargo:rerun-if-changed=src/import.rs");
    println!("cargo:rerun-if-changed=src/import/fastnlo.cpp");
    println!("cargo:rerun-if-changed=src/import/fastnlo.hpp");
}

#[cfg(not(feature = "fastnlo"))]
fn main() {}
