# Contribution guide

## Rust

- Before you commit, first make sure that the whole project including
  libraries, programs, documentation, and tests compile. Use
    - `cargo doc --no-deps --open`, to generate and check the documentation,
    - `cargo test`, to make sure all tests and doctests are still passing,
    - `cargo clippy`, to run the linter and to fix important warnings,
    - and finally `cargo fmt`, to format your code properly.
- Make sure to keep `CHANGELOG.md` up-to-date.
- Make sure not to use Rust features newer than the specified minimum supported
  Rust Version (MSRV), which is documented in the [README](README.md)

## Git

- When you commit, make sure the commit message is written properly. This
  blogpost explains it nicely: <https://chris.beams.io/posts/git-commit/>.
- Whenever possible, prefer rebase over merge.

## Making a new release

Run

    ./make_release 0.5.4

and replace `0.5.4` with a version string, *not* including `v` at the start.
The version strings must adhere to [Semantic
Versioning](https://semver.org/spec/v2.0.0.html).

This will take care of almost everything: the C, Python and Rust interfaces and
their documentation. After some time also a new [Conda
package](https://github.com/conda-forge/pineappl-feedstock) will be generated,
for which the pull request will have to be accepted manually though.
