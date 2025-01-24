#!/bin/bash

set -eou pipefail

./download-test-data.sh

cd ..

# we compile with different flags and don't want to destroy the other target directory
export CARGO_TARGET_DIR="$(mktemp -d)"

# relative paths sometimes don't work, so use an absolute path
dir="${CARGO_TARGET_DIR}"/debug/doctestbins

# `-C link-dead-code` is needed to prevent 'warning: XX functions have mismatched data' warnings
export RUSTFLAGS="-Cinstrument-coverage -Clink-dead-code"
export RUSTDOCFLAGS="-Cinstrument-coverage -Z unstable-options --persist-doctests ${dir}"

# -Z doctest-in-workspace is enabled by default starting from 1.72.0
cargo test -Z doctest-in-workspace --features=applgrid,evolve,fastnlo,fktable 2> >(tee stderr 1>&2)
# from https://stackoverflow.com/a/51141872/812178
sed -i 's/\x1B\[[0-9;]\{1,\}[A-Za-z]//g' stderr

find . -name '*.profraw' -exec $(rustc --print target-libdir)/../bin/llvm-profdata merge -sparse -o pineappl.profdata {} +
( sed -nE 's/[[:space:]]+Running( unittests|) [^[:space:]]+ \(([^)]+)\)/\2/p' stderr && echo "${dir}"/*/rust_out | tr ' ' "\n" ) | \
    xargs printf ' --object %s' | \
    xargs $(rustc --print target-libdir)/../bin/llvm-cov show \
        --ignore-filename-regex='/.cargo/registry' \
        --ignore-filename-regex='rustc' \
        --ignore-filename-regex='pineappl/tests' \
        --ignore-filename-regex='pineappl_capi' \
        --ignore-filename-regex='pineappl_cli/tests' \
        --instr-profile=pineappl.profdata \
        --object "${CARGO_TARGET_DIR}"/debug/pineappl \
        --format html \
        --output-dir cov \
        -Xdemangler=rustfilt

# remove merged profile data and standard error output
rm pineappl.profdata stderr
# remove profile data
find . -name '*.profraw' -delete
# remove build directory
rm -rf "${CARGO_TARGET_DIR}"
