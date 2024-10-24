# Contribution guide

## Rust

- Before you commit, make sure that you have [pre-commit](https://pre-commit.com/)
  installed. This will ensure that the code is formatted correctly and that
  it compiles properly. Also, check if your changes introduce any new linter
  warnings by running `cargo clippy`.
- Make sure to keep `CHANGELOG.md` up-to-date.
- Make sure not to use Rust features newer than the specified minimum supported
  Rust Version (MSRV), which is documented in the [README](README.md). You can
  use `cargo-msrv` to check the crates. However, the Github CI also checks this.
- Make sure to follow the [Rust API Guidelines]

[Rust API Guidelines]: https://rust-lang.github.io/api-guidelines/checklist.html

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

### Writing tests that need test data

- if you write a test that needs test data (grids, EKOs, etc.) store them at
  <https://data.nnpdf.science/pineappl/test-data/>. Ask one of the maintainers
  to upload the data for you if you don't have access to this location. Then
  add a line to `maintainer/generate-coverage.sh` that downloads the data with
  `wget` and a similar line to `.github/workflows/rust.yml` that downloads the
  data with `curl`. To make Github refresh the cached test data when running
  the CI, increase the integer `XX` in the line `key: test-data-vXX` by one.

## Git

- When you commit, make sure the commit message is written properly. This
  blogpost explains it nicely: <https://chris.beams.io/posts/git-commit/>.
- Whenever you have unpushed local commits that are behind `origin/master`, use
  `git pull --rebase` to rebase them
- When editing Github workflow files, use a separate branch, because usually
  many commits are needed to get something working. When merging this branch
  into `master` (or any other branch), squash-merge the commits; the exact
  history in this case is not important

## Making a new release

In the `maintainers` directory run

    ./make_release 0.5.4

and replace `0.5.4` with a version string, *not* including `v` at the start.
The version strings must adhere to [Semantic Versioning].

This will take care of almost everything: the C, Python and Rust interfaces and
their documentation. After some time also a new [Conda package] will be
generated, for which the pull request will have to be accepted manually though.

[Semantic Versioning]: https://semver.org/spec/v2.0.0.html
[Conda package]: https://github.com/conda-forge/pineappl-feedstock

## Updating the CI's container

To update the software the CI runs with, modify the files in
`maintainer/pineappl-ci`. See `maintainer/README.md` for a description of what
these files do. To generate a new container, you need to manually run the
[Container action] from the branch in which you modified the container files.
After the container has been generated, all following commits in *every* branch
will use the new container.

[Container action]: https://github.com/NNPDF/pineappl/actions/workflows/container.yml
