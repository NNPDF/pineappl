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

curl -s -LJ "${base_url}/v${version}/pineappl_cli-${target}.tar.gz" \
    | tar xzf - -C "${prefix}"

if command -v pineappl >/dev/null; then
    path="$(command -v pineappl)"

    if [ "${path}" != "${prefix}"/bin/pineappl ]; then
        echo
        echo "Warning: Your PATH evironment variable isn't properly set."
        echo "It appears a different installation of PineAPPL is found:"
        echo
        echo "  ${path}"
        echo
        echo "Remove this installation or reorder your PATH"
    fi
else
    echo
    echo "Warning: Your PATH environment variable isn't properly set."
    echo "Try adding"
    echo
    echo "  export PATH=${prefix}\"/bin:\${PATH}\""
    echo
    echo "to your shell configuration file"
fi
