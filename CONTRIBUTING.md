# Contribution guide

- Before you commit, make sure to first run the following commands and fix all
  warnings and errors, in each step, before running the next command:
    - `cargo build`, to build the project,
    - `cargo clippy`, to see whether the linter has any good suggestions,
    - `cargo test`, to make sure all tests and doctests are still passing,
    - `cargo fmt`, to format your code properly,
    - `cargo doc --open`, to generate the documentation. Make sure the
      documentation is formatted properly.
