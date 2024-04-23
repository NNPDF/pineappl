#!/bin/sh

# WARNING: do not commit changes to this file unless you've checked it against
# `shellcheck` (https://www.shellcheck.net/); run `shellcheck install-capi.sh`
# to make sure this script is POSIX shell compatible; we cannot rely on bash
# being present

set -eu

prefix=
version=

while [ $# -gt 0 ]; do
    case $1 in
        --version)
            version=$2
            shift
            shift
            ;;
        --version=*)
            version=${1#--version=}
            shift
            ;;
        --prefix)
            prefix=$2
            shift
            shift
            ;;
        --prefix=*)
            prefix=${1#--prefix=}
            shift
            ;;
        --target)
            target=$2
            shift
            shift
            ;;
        --target=*)
            target=${1#--target=}
            shift
            ;;
        *)
            echo "Error: argument '$1' unknown"
            exit 1
            ;;
    esac
done

if [ -z ${target+x} ]; then
    case $(uname -m):$(uname -s) in
        arm64:Darwin)
            target=aarch64-apple-darwin;;
        x86_64:Darwin)
            target=x86_64-apple-darwin;;
        x86_64:Linux)
            target=x86_64-unknown-linux-gnu;;
        *)
            echo "Error: unknown target, uname = '$(uname -a)'"
            exit 1;;
    esac
fi

# if no prefix is given, prompt for one
if [ -z "${prefix}" ]; then
    # read from stdin (`<&1`), even if piped into a shell
    printf "Enter installation path: "
    read -r <&1 prefix
    echo
fi

# we need the absolute path; use `eval` to expand possible tilde `~`
eval mkdir -p "${prefix}"
eval cd "${prefix}"
prefix=$(pwd)
cd - >/dev/null

# if no version is given, use the latest version
if [ -z "${version}" ]; then
    version=$(curl -s https://api.github.com/repos/NNPDF/pineappl/releases/latest | \
        sed -n 's/[ ]*"tag_name"[ ]*:[ ]*"v\([^"]*\)"[ ]*,[ ]*$/\1/p')
fi

base_url=https://github.com/NNPDF/pineappl/releases/download

echo "prefix:  ${prefix}"
echo "target:  ${target}"
echo "version: ${version}"

curl -s -LJ "${base_url}/v${version}/pineappl_capi-${target}.tar.gz" \
    | tar xzf - -C "${prefix}"

# instead of `sed` and `mv` we could use `sed -i`, but on Mac it doesn't work as expected from GNU sed
sed "s:prefix=/:prefix=${prefix}:" "${prefix}"/lib/pkgconfig/pineappl_capi.pc > \
    "${prefix}"/lib/pkgconfig/pineappl_capi.pc.new
mv "${prefix}"/lib/pkgconfig/pineappl_capi.pc.new "${prefix}"/lib/pkgconfig/pineappl_capi.pc

pcbin=

if command -v pkg-config >/dev/null; then
    pcbin=$(command -v pkg-config)
elif command -v pkgconf >/dev/null; then
    pcbin=$(command -v pkgconf)
else
    echo
    echo "Warning: neither \`pkg-config\` nor \`pkgconf\` not found. At least one is needed for the CAPI to be found"
    exit 1
fi

# check whether the library can be found
if "${pcbin}" --libs pineappl_capi >/dev/null 2>/dev/null; then
    prefix_lib=$(cd "${prefix}"/lib && pwd)
    found_lib=$("${pcbin}" --keep-system-libs --libs-only-L pineappl_capi | sed 's/-L[[:space:]]*//')

    if [ "${prefix_lib}" != "${found_lib}" ]; then
        echo
        echo "Warning: Your PKG_CONFIG_PATH environment variable isn't properly set."
        echo "It appears a different installation of PineAPPL is found:"
        echo
        echo "  ${found_lib}"
        echo
        echo "Remove this installation or reorder your PKG_CONFIG_PATH"
    fi
else
    echo
    echo "Warning: Your PKG_CONFIG_PATH environment variable isn't properly set."
    echo "Try adding"
    echo
    echo "  export PKG_CONFIG_PATH=${prefix}/lib/pkgconfig${PKG_CONFIG_PATH:+:\"\${PKG_CONFIG_PATH\}\"}"
    echo
    echo "to your shell configuration file"
fi
