#!/bin/bash

set -euo pipefail

# install rustup
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

for version in ${RUST_V}; do
    # the last command will be the default
    rustup default ${version}
    # install LLVM tools needed for code coverage
    rustup component add llvm-tools-preview
done

# # install Fortran compiler
# apt update
# apt install gfortran -y

# install cargo-c needed for the CAPI
cargo install cargo-c --version ${CARGOC_V} --features=vendored-openssl

# remove files generated by cargo
rm -r /usr/local/cargo/registry

# install LHAPDF
curl "https://lhapdf.hepforge.org/downloads/?f=LHAPDF-${LHAPDF_V}.tar.gz" | tar xzf -
cd LHAPDF-${LHAPDF_V}
./configure --disable-python --disable-static
make -j
make install
ldconfig

cd ..

# install PDF sets
for pdf in NNPDF31_nlo_as_0118_luxqed NNPDF40_nnlo_as_01180 NNPDF40_nlo_as_01180; do
    curl "https://lhapdfsets.web.cern.ch/current/${pdf}.tar.gz" | tar xzf - -C /usr/local/share/LHAPDF
done

# install APPLgrid
curl "https://applgrid.hepforge.org/downloads?f=applgrid-${APPLGRID_V}.tgz" | tar xzf -
cd applgrid-${APPLGRID_V}
./configure --disable-static --without-root
make -j
make install
ldconfig
mkdir -p ${APPL_IGRID_DIR}
cp src/*.h ${APPL_IGRID_DIR}

cd ..

# install fastNLO
curl "https://fastnlo.hepforge.org/code/v25/fastnlo_toolkit-${FASTNLO_V}.tar.gz" | tar xzf -
cd fastnlo_toolkit-${FASTNLO_V}
./configure --disable-static --prefix=/usr/local/
make -j
make install
ldconfig

cd ..
