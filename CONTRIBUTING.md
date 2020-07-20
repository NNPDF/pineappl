# Contribution guide

- Before you commit, first make sure that the whole project including
  libraries, programs, documentation, and tests compile. Use
    - `cargo doc --no-deps --open`, to generate and check the documentation,
    - `cargo test`, to make sure all tests and doctests are still passing,
    - `cargo clippy`, to run the linter and to fix important warnings,
    - and finally `cargo fmt`, to format your code properly.
- Make sure to keep `CHANGELOG.md` up-to-date.
- When you commit, make sure the commit message is written properly. This
  blogpost explains it nicely: https://chris.beams.io/posts/git-commit/.
