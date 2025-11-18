fn main() {

    let out = std::env::var("OUT_DIR").unwrap();

    cc::Build::new()
        .file("src/fx2.c")
        .flag("-O2")              
        .flag("-ffp-contract=off") 
        .compile("fx2");

    println!("cargo:rustc-link-lib=static=fx2");
    println!("cargo:rustc-link-search=native={}", out);
    println!("cargo:rerun-if-changed=src/fx2.c");
}
