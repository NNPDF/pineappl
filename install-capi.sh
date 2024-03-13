#!/bin/sh

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
    if command -v gcc >/dev/null; then
        target=$(gcc -dumpmachine)

        case ${target} in
            arm64-apple-darwin*) target=aarch64-apple-darwin;;
            x86_64-*-linux-gnu | x86_64-linux-gnu | x86_64-*-linux) target=x86_64-unknown-linux-gnu;;
            x86_64-apple-darwin*) target=x86_64-apple-darwin;;
            *) echo "Error: target '${target}' unknown."; exit 1;;
        esac
    else
        echo "Error: target unknown, try setting it manually with '--target'"
        exit 1
    fi
fi

# if no prefix is given, prompt for one
if [ -z ${prefix} ]; then
    # read from stdin (`<&1`), even if piped into a shell
    read -p "Enter installation path: " <&1 prefix
fi

if [ ! -d "${prefix}" ]; then
    mkdir -p "${prefix}"
fi

# if no version is given, use the latest version
if [ -z ${version} ]; then
    version=$(curl -s https://api.github.com/repos/NNPDF/pineappl/releases/latest | \
        sed -n 's/[ ]*"tag_name"[ ]*:[ ]*"v\([^"]*\)"[ ]*,[ ]*$/\1/p')
fi

base_url=https://github.com/NNPDF/pineappl/releases/download

echo "prefix:  ${prefix}"
echo "target:  ${target}"
echo "version: ${version}"

# we need the absolute path
cd "${prefix}"
prefix=$(pwd)
cd - >/dev/null

curl -s -LJ "${base_url}"/v${version}/pineappl_capi-${target}.tar.gz \
    | tar xzf - -C "${prefix}"

# instead of `sed` and `mv` we could use `sed -i`, but on Mac it doesn't work as expected from GNU sed
sed s:prefix=/:prefix=${prefix}: "${prefix}"/lib/pkgconfig/pineappl_capi.pc > \
    "${prefix}"/lib/pkgconfig/pineappl_capi.pc.new
mv "${prefix}"/lib/pkgconfig/pineappl_capi.pc.new "${prefix}"/lib/pkgconfig/pineappl_capi.pc

if command -v pkg-config >/dev/null; then
    if ! pkg-config --libs pineappl_capi >/dev/null; then
        echo
        echo "Warning: Your PKG_CONFIG_PATH environment variable isn't properly set. Try adding"
        echo "  export PKG_CONFIG_PATH=${prefix}/lib/pkgconfig"
        echo "to your shell configuration file"
    fi
else
    echo
    echo "Warning: \`pkg-config\` binary not found. Without it the CAPI may not be found."
fi
