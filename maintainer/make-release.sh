#!/bin/bash

set -euo pipefail

crates=(
    pineappl
    pineappl_applgrid
    pineappl_fastnlo
    pineappl_capi
    pineappl_cli
    pineappl_py
    xtask
)

features=(
    applgrid
    evolve
    fastnlo
    fktable
)

main=master
this_branch=$(git rev-parse --abbrev-ref HEAD)

cd ..

if [[ $# != 1 ]]; then
    echo "No version number given."
    exit 1
fi

version=$1

prerelease=$(echo ${version} | perl -pe 's/^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)(?:-((?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?(?:\+([0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$/\4/')

if [[ $(echo ${version} | grep -oP '^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)(?:-((?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?(?:\+([0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$') != ${version} ]]; then
    echo "Version string incorrect."
    exit 1
fi

# in branches that are not master we only allow prereleases
if [[ ${this_branch} != ${main} ]] && [[ ${prerelease} == "" ]]; then
    echo "Ordinary releases are only allowed in the '${main}' branch."
    echo "If you really want to make a release from '${this_branch}', consider making a prerelease."
    exit 1
fi

for crate in ${crates[@]}; do
    if [[ -n $(git status ${crate} --porcelain) ]]; then
        echo "This repository isn't clean. Make sure to add or delete the corresponding files."
        exit 1
    fi
done

echo ">>> Updating version strings ..."

# we don't want to create a changelog entry for prereleases, which are solely
# for internal testing purposes
if [[ ${prerelease} == "" ]]; then
    sed -i \
        -e "s:\(## \[Unreleased\]\):\1\n\n## [${version}] - $(date +%d/%m/%Y):" \
        -e "s:\[Unreleased\]\(\: https\://github.com/NNPDF/pineappl/compare/v\)\(.*\)...HEAD:[Unreleased]\1${version}...HEAD\n[${version}]\1\2...v${version}:" \
        CHANGELOG.md
    git add CHANGELOG.md
fi

# the '.' is needed because we also need to modify the workspace
for crate in . ${crates[@]}; do
    sed -i \
        -e "s:^version = \".*\":version = \"${version}\":" \
        -e "s:^\(pineappl = .*\)version = \".*\":\1version = \"=${version}\":" \
        -e "s:^\(pineappl_applgrid = .*\)version = \".*\":\1version = \"=${version}\":" \
        -e "s:^\(pineappl_cli = .*\)version = \".*\":\1version = \"=${version}\":" \
        -e "s:^\(pineappl_fastnlo = .*\)version = \".*\":\1version = \"=${version}\":" \
        ${crate}/Cargo.toml
    git add ${crate}/Cargo.toml
done

echo ">>> Updating Cargo.lock ..."

for crate in "${crates[@]}"; do cargo pkgid path+file://$(cd ../"${crate}" && pwd); done | xargs printf " -p %s" | xargs cargo update
git add Cargo.lock

echo ">>> Testing if 'pineappl' can be published ..."

cd pineappl
cargo publish --dry-run
cd ..

echo ">>> Commiting and pushing changes ..."

git commit -m "Release v${version}"
git tag -a v${version} -m v${version}
git push --follow-tags
