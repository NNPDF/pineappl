use super::helpers::Subcommand;
use anyhow::{anyhow, Result};
use clap::{Parser, ValueHint};
use std::path::PathBuf;

// the indirection of two modules instead of one is needed to hide the C++ bridge behind the
// feature gate
#[cfg(feature = "import-fastnlo")]
mod fastnlo {
    #[cxx::bridge]
    pub mod ffi {
        unsafe extern "C++" {
            include!("pineappl_cli/src/import/fastnlo.hpp");

            unsafe fn this_would_be_main(input: *const c_char, output: *const c_char) -> i32;
        }
    }
}

/// Converts fastNLO tables to PineAPPL grids.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the converted grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    output: PathBuf,
}

impl Subcommand for Opts {
    #[cfg(feature = "import-fastnlo")]
    fn run(&self) -> Result<()> {
        use std::ffi::CString;
        use std::os::unix::ffi::OsStrExt;

        // TODO: this is a correct way of converting strings, but only works on UNIX
        let input = CString::new(self.input.as_path().as_os_str().as_bytes()).unwrap();
        let output = CString::new(self.output.as_path().as_os_str().as_bytes()).unwrap();

        let errcode = unsafe { fastnlo::ffi::this_would_be_main(input.as_ptr(), output.as_ptr()) };

        if errcode == 0 {
            Ok(())
        } else {
            Err(anyhow!("something went wrong"))
        }
    }

    #[cfg(not(feature = "import-fastnlo"))]
    fn run(&self) -> Result<()> {
        Err(anyhow!(
            "you need to install `pineappl` with feature `import-fastnlo`"
        ))
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["import", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }
}
