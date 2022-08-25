#!/bin/bash

set -e

version=$1

if [[ -z ${version} ]]; then
    echo "No version number given."
    exit 1
fi

if [[ ! ${version} =~ ^[0-9]+.[0-9]+.[0-9]+$ ]]; then
    echo "Version string incorrect."
    exit 1
fi

if ! which gh >/dev/null; then
    echo "Didn't find the \`gh\` binary."
    exit 1
fi

if ! gh auth status 2>/dev/null; then
    echo "Couldn't connect to the github repository."
    exit 1
fi

if ! cargo msrv --help >/dev/null; then
    echo "Didn't find \`msrv\` applet of \`cargo\`."
    exit 1
fi

if ! cargo msrv --min 1.56.1 --max 1.56.1 >/dev/null; then
    echo "Minimum supported Rust version doesn't match avertised one."
    exit
fi

#if [[ -n $(git status --porcelain) ]]; then
#    echo "This repository isn't clean. Make sure to add or delete the corresponding files."
#    exit 1
#fi

#if [[ ]]; then
#    echo "You're not on master."
#    exit 1
#fi

echo ">>> Testing release configuration with default features ..."

cargo build --release
cargo test --release

echo ">>> Testing release configuration with \`applgrid\` feature ..."

cargo build --release --features=applgrid
cargo test --release --features=applgrid

echo ">>> Testing release configuration with \`fastnlo\` feature ..."

cargo build --release --features=fastnlo
cargo test --release --features=fastnlo

echo ">>> Testing release configuration with \`fktable\` feature ..."

cargo build --release --features=fktable
cargo test --release --features=fktable

echo ">>> Testing if 'pineappl' can be published ..."

cd pineappl
cargo publish --dry-run
cd ..

echo ">>> Updating version strings ..."

sed -i \
    -e "s:\(## \[Unreleased\]\):\1\n\n## [${version}] - $(date +%d/%m/%Y):" \
    -e "s:\[Unreleased\]\(\: https\://github.com/N3PDF/pineappl/compare/v\)\(.*\)...HEAD:[Unreleased]\1${version}...HEAD\n[${version}]\1\2...v${version}:" \
    CHANGELOG.md

sed -i \
    -e "s:^version = \".*\":version = \"${version}\":" \
    -e "s:^\(pineappl = .*\)version = \".*\":\1version = \"${version}\":" \
    pineappl{,_capi,_cli,_py}/Cargo.toml

echo ">>> Commiting and pushing changes ..."

git commit -a -m "Release v${version}"
git tag -a v${version} -m v${version}
git push --follow-tags

echo ">>> Publishing crate 'pineappl' ..."

cd pineappl
cargo publish

echo "Waiting the 'pineappl' crate to become available on crates.io ..."

sleep 60

echo ">>> Publishing crate 'pineappl_capi' ..."

cd ../pineappl_capi
cargo publish

echo ">>> Publishing crate 'pineappl_cli' ..."

cd ../pineappl_cli
cargo publish
cd ..

echo ">>> Making a release on github"

# extract the previous version number
old_version=$(sed -n 's/^## \[\(.*\)\] - .*/\1/p' CHANGELOG.md | tail +2 | head -n 1)

# extract news for the current version from the changelog file, dismissing
# empty lines at the start and the end
news=$(sed -n "/\\[${version}\\]/, /\\[${old_version}\\]/{ /\\[${version}\\]/! { /\\[${old_version}\\]/! p } }" \
    CHANGELOG.md | sed -e :a -e '/./,$!d;/^\n*$/{$d;N;};/\n$/ba')

gh release create v${version} -n "${news}"
