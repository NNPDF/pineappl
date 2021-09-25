# Contribution guide

## Rust

- Before you commit, first make sure that the whole project including
  libraries, programs, documentation, and tests compile. Use
    - `cargo doc --no-deps --open`, to generate and check the documentation,
    - `cargo test`, to make sure all tests and doctests are still passing,
    - `cargo clippy`, to run the linter and to fix important warnings,
    - and finally `cargo fmt`, to format your code properly.
- Make sure to keep `CHANGELOG.md` up-to-date.

## Git

- When you commit, make sure the commit message is written properly. This
  blogpost explains it nicely: https://chris.beams.io/posts/git-commit/.
- Whenever possible, prefer rebase over merge

## Releasing a tagged version

Please follow these steps to make a release:

1) Make sure to correctly set the version number you want to release in
    - `pineappl/Cargo.toml`,
    - `pineappl_capi/Cargo.toml`,
    - `pineappl_cli/Cargo.toml` and
    - `pineappl_py/Cargo.toml`,
   both in the `version` field and the dependency on `pineappl`, if that exists.
2) Update `CHANGELOG.md` by adding the version you are about to release before
   the 'Unreleased' header and update the necessary links.
3) Commit the changes with the message 'Release vX.X.X'
4) Tag the release, and push both the commit and the tag to github.com
5) Create a release on github.com, using the entries of `CHANGELOG.md` for this
   version
6) Immediately create another release increasing the version number for all
   subsequent commits; push it immediately
