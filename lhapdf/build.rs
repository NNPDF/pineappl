fn main() {
    let lhapdf = pkg_config::Config::new()
        .atleast_version("6")
        .probe("lhapdf")
        .unwrap();
    let mut cpp_build = cpp_build::Config::new();

    for include_path in lhapdf.include_paths {
        cpp_build.include(include_path);
    }

    cpp_build.build("src/lib.rs");

    for lib_path in lhapdf.link_paths {
        println!("cargo:rustc-link-search={}", lib_path.to_str().unwrap());
    }

    for lib in lhapdf.libs {
        println!("cargo:rustc-link-lib={}", lib);
    }
}
