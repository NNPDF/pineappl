# Contribution guide

## Rust

- Before you commit, make sure that your code compiles with `cargo check` and
  that it has been formatted properly; `cargo fmt` does that for you.
- Make sure to keep `CHANGELOG.md` up-to-date.
- Make sure not to use Rust features newer than the specified minimum supported
  Rust Version (MSRV), which is documented in the [README](README.md). You can
  use `cargo-msrv` to check the crates. However, the Github CI also checks this.
- Make sure to follow the [Rust API
  Guidelines](https://rust-lang.github.io/api-guidelines/checklist.html)

### Increasing the minimum supported Rust version (MSRV)

Do not change the MSRV for releases with increased patch version number. When
increasing the MSRV make sure to set it everywhere to the same value:

1. first update the `pineappl-ci` container, by
   - adding the new MSRV to the variable `RUST_V`,
   - making sure the nightly version is the last entry and
   - leaving in the previous MSRV to not break the CI in between the transition
     from it
2. commit the previous changes and manually run the `Container` Github action
3. next, update the MSRV in the following files:
   - in the top-level `Cargo.toml`; all other projects in the workspace should
     inherit the setting in their respective `Cargo.toml` files
   - in `README.md` and `docs/installation.md`
   - in all Github workflows (`.github/workflows/`)
   - in `.readthedocs.yml` update the value of the `rust` field and make sure
     [RTD supports it](https://docs.readthedocs.io/en/stable/config-file/v2.html#build-tools-rust)
   - in `make_release.sh` update the `cargo msrv` call
4. commit the previous changes and push them *after* the container created by
   step 2 is ready

### Coding guidelines

- avoid the use of indices whenever possible; use `Iterator` instead.
- use the `unwrap` methods whenever a panic would signal a bug in the program,
  and use `Result` instead if errors should be propagated down to the user.
  When using `unwrap`, document the nature of the bug if a panic happens with a
  comment of the form: `// UNWRAP: ...`.
- in APIs prefer `unwrap_or_else(|| unreachable!())` over `unwrap` whenever
  this avoids the clippy warning that a Panic section is missing. Also document
  this with `// UNWRAP: ...`

## Git

- When you commit, make sure the commit message is written properly. This
  blogpost explains it nicely: <https://chris.beams.io/posts/git-commit/>.
- Whenever possible, prefer rebase over merge.

## Making a new release

In the `maintainers` directory run

    ./make_release 0.5.4

and replace `0.5.4` with a version string, *not* including `v` at the start.
The version strings must adhere to [Semantic
Versioning](https://semver.org/spec/v2.0.0.html).

This will take care of almost everything: the C, Python and Rust interfaces and
their documentation. After some time also a new [Conda
package](https://github.com/conda-forge/pineappl-feedstock) will be generated,
for which the pull request will have to be accepted manually though.
