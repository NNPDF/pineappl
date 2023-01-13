# Contribution guide

## Rust

- Before you commit, make sure that your code compiles with `cargo check` and
  that it has been formatted properly; `cargo fmt` does that for you.
- Make sure to keep `CHANGELOG.md` up-to-date.
- Make sure not to use Rust features newer than the specified minimum supported
  Rust Version (MSRV), which is documented in the [README](README.md). You can
  use `cargo-msrv` to check the crates. However, the Github CI also checks this.

### Coding guidelines

- avoid the use of indices whenever possible; use `Iterator` instead.
- use the `unwrap` methods whenever a panic would signal a bug in the program,
  and use `Result` instead if errors should be propagated down to the user.
- in APIs prefer `unwrap_or_else(|| unreachable!())` over `unwrap` whenever
  this avoids the clippy warning that a Panic section is missing

## Git

- When you commit, make sure the commit message is written properly. This
  blogpost explains it nicely: <https://chris.beams.io/posts/git-commit/>.
- Whenever possible, prefer rebase over merge.

## Increasing the minimum supported Rust version (MSRV)

Do not change the MSRV for releases with increased patch version number. When
increasing the MSRV make sure to set it everywhere to the same value:

- update `rust-version` in each crate of the workspace
- update the MSRV in `README.md` and `docs/installation.md`
- update the MSRV in all Github workflows (`.github/workflows/`)
- update `rust` in `.readthedocs.yml` and make sure RTD supports it

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
