#!/bin/bash

set -eou pipefail

download() {
    if command -v wget >/dev/null; then
        wget --no-verbose --no-clobber -P test-data "$1"
    elif command -v curl >/dev/null; then
        test -d test-data || mkdir test-data
        cd test-data
        curl -s -C - -O "$1"
        cd ..
    else
        echo "neither 'wget' nor 'curl' found, exiting."
        exit 1
    fi
}

files=(
    'https://data.nnpdf.science/dy_high_mass/CMS_DY_14TEV_MLL_6000_COSTH.pineappl.lz4'
    'https://data.nnpdf.science/dy_high_mass/NNPDF_DY_14TEV_BSM_AFB.pineappl.lz4'
    'https://data.nnpdf.science/pineappl/test-data/ATLASWPT11-Wplus_tot.appl'
    'https://data.nnpdf.science/pineappl/test-data/CMS_TTB_8TEV_2D_TTM_TRAP_TOT-opt.pineappl.lz4'
    'https://data.nnpdf.science/pineappl/test-data/CMS_TTB_8TEV_2D_TTM_TRAP_TOT.tar'
    'https://data.nnpdf.science/pineappl/test-data/E906nlo_bin_00.pineappl.lz4'
    'https://data.nnpdf.science/pineappl/test-data/E906nlo_bin_00.tar'
    'https://data.nnpdf.science/pineappl/test-data/FKTABLE_CMSTTBARTOT8TEV-TOPDIFF8TEVTOT.pineappl.lz4'
    'https://data.nnpdf.science/pineappl/test-data/FKTABLE_STAR_WMWP_510GEV_WM-AL-POL.pineappl.lz4'
    'https://data.nnpdf.science/pineappl/test-data/FK_ATLASTTBARTOT13TEV.dat'
    'https://data.nnpdf.science/pineappl/test-data/FK_POSXDQ.dat'
    'https://data.nnpdf.science/pineappl/test-data/GRID_DYE906R_D_bin_1.pineappl.lz4'
    'https://data.nnpdf.science/pineappl/test-data/GRID_STAR_WMWP_510GEV_WP-AL-POL.pineappl.lz4'
    'https://data.nnpdf.science/pineappl/test-data/LHC8-Mtt-HT4-173_3-bin1.tab.gz'
    'https://data.nnpdf.science/pineappl/test-data/LHCBWZMU7TEV_PI_part1.appl'
    'https://data.nnpdf.science/pineappl/test-data/LHCB_DY_8TEV.pineappl.lz4'
    'https://data.nnpdf.science/pineappl/test-data/LHCB_DY_8TEV.tar'
    'https://data.nnpdf.science/pineappl/test-data/LHCB_WP_7TEV.tar'
    'https://data.nnpdf.science/pineappl/test-data/EKO_LHCB_WP_7TEV.txt'
    'https://data.nnpdf.science/pineappl/test-data/LHCB_WP_7TEV_old.pineappl.lz4'
    'https://data.nnpdf.science/pineappl/test-data/LHCB_WP_7TEV_opt.pineappl.lz4'
    'https://data.nnpdf.science/pineappl/test-data/LHCB_WP_7TEV_v2.tar'
    'https://data.nnpdf.science/pineappl/test-data/LHCB_WP_7TEV_v2_xif_2.tar'
    'https://data.nnpdf.science/pineappl/test-data/NJetEvents_0-0-2.tab.gz'
    'https://data.nnpdf.science/pineappl/test-data/NNPDF_POS_F2D_40.pineappl.lz4'
    'https://data.nnpdf.science/pineappl/test-data/NUTEV_CC_NU_FE_SIGMARED.pineappl.lz4'
    'https://data.nnpdf.science/pineappl/test-data/NUTEV_CC_NU_FE_SIGMARED.tar'
    'https://data.nnpdf.science/pineappl/test-data/SIHP-PP-POLARIZED-STAR-NLO.pineappl.lz4'
    'https://data.nnpdf.science/pineappl/test-data/STAR_WMWP_510GEV_WM-AL-POL.pineappl.lz4'
    'https://data.nnpdf.science/pineappl/test-data/STAR_WMWP_510GEV_WM-AL-POL_PolPDF.tar'
    'https://data.nnpdf.science/pineappl/test-data/STAR_WMWP_510GEV_WM-AL-POL_UnpolPDF.tar'
    'https://data.nnpdf.science/pineappl/test-data/ZEUS_2JET_319GEV_374PB-1_DIF_ETQ2_BIN6.pineappl.lz4'
    'https://data.nnpdf.science/pineappl/test-data/ZEUS_2JET_319GEV_374PB-1_DIF_ETQ2_BIN6.tar'
    'https://data.nnpdf.science/pineappl/test-data/LHCB_WP_8TEV.pineappl.lz4'
    'https://data.nnpdf.science/pineappl/test-data/LHCB_WP_8TEV.tar'
    'https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-atlas-dijets-fnlo-arxiv-1312.3524/grids/applfast-atlas-dijets-fnlo-arxiv-1312.3524-xsec000.tab.gz'
    'https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-dijets-appl-arxiv-0010054/grids/applfast-h1-dijets-appl-arxiv-0010054-xsec000.appl'
    'https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-incjets-fnlo-arxiv-0706.3722/grids/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz'
    'https://ploughshare.web.cern.ch/ploughshare/db/atlas/atlas-atlas-wpm-arxiv-1109.5141/grids/atlas-atlas-wpm-arxiv-1109.5141-xsec001.appl'
)

cd ..

for file in "${files[@]}"; do
    download "${file}"
done
