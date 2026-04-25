#!/bin/bash

APPLGRID_V=1.6.36
FASTNLO_V=2.6
FASTNLO_R=69d87cf4
LHAPDF_V=6.5.4
ZLIB_V=1.3.1

# if we're in a container, this variable is already set by the container
if [[ ! -v APPL_IGRID_DIR ]]; then
    export APPL_IGRID_DIR=/usr/local/src/applgrid/src
fi

# downloading files is fragile, so fail early
urls=(
    "https://lhapdf.hepforge.org/downloads/?f=LHAPDF-${LHAPDF_V}.tar.gz"
    "https://www.zlib.net/fossils/zlib-${ZLIB_V}.tar.gz"
    "https://applgrid.hepforge.org/downloads/applgrid-${APPLGRID_V}.tgz"
    # "https://fastnlo.hepforge.org/code/v25/fastnlo_toolkit-${FASTNLO_V}.tar.gz"
)

# download fastNLO from Git to support PineAPPL -> fastNLO export
git clone https://gitlab.etp.kit.edu/qcd-public/fastNLO.git

for url in "${urls[@]}"; do
    curl -fsSL "${url}" | tar xzf -
done

# install LHAPDF
cd "LHAPDF-${LHAPDF_V}"
if [[ ${static} == yes ]]; then
    # compile static libraries with PIC to make statically linking PineAPPL's CLI work
    # see also https://users.rust-lang.org/t/why-does-crelocation-model-dynamic-no-pic-help-although-it-shouldnt/109012
    ./configure --disable-python --disable-shared --with-pic=yes
else
    ./configure --disable-python
fi
make -j V=1
make install
ldconfig
cd ..

# install zlib
cd "zlib-${ZLIB_V}"
if [[ ${static} == yes ]]; then
    CFLAGS=-fPIC ./configure --static --prefix=/usr/local
else
    ./configure --prefix=/usr/local
fi
make -j
make install
ldconfig
cd ..

# install APPLgrid
cd "applgrid-${APPLGRID_V}"
patch -l -p0 <<EOF
--- src/combine.cxx	2024-04-23 16:35:27.000000000 +0200
+++ src/combine.cxx.new	2024-07-06 12:29:12.813303074 +0200
@@ -56,12 +56,6 @@
 }


-double integral( appl::TH1D* h ) {
-  double d = 0;
-  for ( int i=0 ; i<h->GetNbinsX() ; i++ ) d += h->GetBinContent(i+1);
-  return d;
-}
-

 void print( appl::TH1D* h ) {
   for ( int i=1 ; i<=h->GetNbinsX() ; i++ ) std::cout << h->GetBinContent(i) << " ";
EOF
if [[ ${static} == yes ]]; then
    # compile static libraries with PIC to make statically linking PineAPPL's CLI work
    ./configure --without-root --disable-shared --with-pic=yes
else
    ./configure --without-root
fi
make -j
make install
ldconfig

# PineAPPL's APPLgrid interface needs some information not in the public headers
mkdir -p "${APPL_IGRID_DIR}"
cp src/*.h "${APPL_IGRID_DIR}"
cd ..

# install fastNLO
cd fastNLO
git checkout "${FASTNLO_R}"
cd "v${FASTNLO_V}/toolkit"
autoreconf -fi
if [[ ${static} == yes ]]; then
    # compile static libraries with PIC to make statically linking PineAPPL's CLI work
    ./configure --prefix=/usr/local/ --disable-shared --with-pic=yes
else
    ./configure --prefix=/usr/local/
fi
make -j V=1
make install
ldconfig
cd ..
