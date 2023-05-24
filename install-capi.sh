#!/bin/sh

set -eu

if command -v gcc >/dev/null; then
    target=$(gcc -dumpmachine)

    case ${target} in
        x86_64-*-linux-gnu | x86_64-linux-gnu | x86_64-*-linux) target=x86_64-unknown-linux-gnu;;
        *) echo "target '${target}' unknown."; exit 1;;
    esac
else
    echo "target unknown."
    exit 1
fi

if [ $# -eq 1 ]; then
    version=$1
else
    # if no version is given, use the latest version
    version=$(curl -s https://api.github.com/repos/NNPDF/pineappl/releases/latest | \
        sed -n 's/[ ]*"tag_name"[ ]*:[ ]*"v\([^"]*\)"[ ]*,[ ]*$/\1/p')
fi

base_url=https://github.com/NNPDF/pineappl/releases/download

echo "target:  ${target}"
echo "version: ${version}"
echo

# read from stdin (`<&1`), even if piped into a shell
read -p "Enter installation path: " <&1 prefix

if test ! -d "${prefix}"; then
    echo "'${prefix}' doesn't exist, can't install into non-existing directories."
    exit 1
fi

# we need the absolute path
cd "${prefix}"
prefix=$(pwd)
cd - >/dev/null

wget --quiet "${base_url}"/v${version}/pineappl_capi-${target}.tar.gz -O- \
    | tar xzf - -C "${prefix}"
sed -i s:prefix=/:prefix=${prefix}: "${prefix}"/lib/pkgconfig/pineappl_capi.pc
