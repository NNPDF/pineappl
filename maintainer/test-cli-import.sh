#!/bin/bash -x

set -euo pipefail

grids=(
    ## Ploughshare

    # Group: applfast
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-atlas-dijets-appl-arxiv-1312.3524/applfast-atlas-dijets-appl-arxiv-1312.3524.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-atlas-dijets-appl-arxiv-1711.02692/applfast-atlas-dijets-appl-arxiv-1711.02692.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-atlas-dijets-fnlo-arxiv-1312.3524/applfast-atlas-dijets-fnlo-arxiv-1312.3524.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-atlas-dijets-fnlo-arxiv-1711.02692/applfast-atlas-dijets-fnlo-arxiv-1711.02692.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-atlas-incjets-appl-arxiv-1410.8857/applfast-atlas-incjets-appl-arxiv-1410.8857.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-atlas-incjets-appl-arxiv-1706.03192/applfast-atlas-incjets-appl-arxiv-1706.03192.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-atlas-incjets-appl-arxiv-1711.02692/applfast-atlas-incjets-appl-arxiv-1711.02692.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-atlas-incjets-fnlo-arxiv-1410.8857/applfast-atlas-incjets-fnlo-arxiv-1410.8857.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-atlas-incjets-fnlo-arxiv-1706.03192/applfast-atlas-incjets-fnlo-arxiv-1706.03192.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-atlas-incjets-fnlo-arxiv-1711.02692/applfast-atlas-incjets-fnlo-arxiv-1711.02692.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-cms-dijets-appl-arxiv-1212.6660/applfast-cms-dijets-appl-arxiv-1212.6660.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-cms-dijets-appl-arxiv-1705.02628/applfast-cms-dijets-appl-arxiv-1705.02628.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-cms-dijets-fnlo-arxiv-1212.6660/applfast-cms-dijets-fnlo-arxiv-1212.6660.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-cms-dijets-fnlo-arxiv-1705.02628/applfast-cms-dijets-fnlo-arxiv-1705.02628.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-cms-incjets-appl-arxiv-1212.6660/applfast-cms-incjets-appl-arxiv-1212.6660.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-cms-incjets-appl-arxiv-1512.06212/applfast-cms-incjets-appl-arxiv-1512.06212.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-cms-incjets-appl-arxiv-1609.05331/applfast-cms-incjets-appl-arxiv-1609.05331.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-cms-incjets-appl-arxiv-2111.10431/applfast-cms-incjets-appl-arxiv-2111.10431.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-cms-incjets-fnlo-arxiv-1212.6660/applfast-cms-incjets-fnlo-arxiv-1212.6660.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-cms-incjets-fnlo-arxiv-1512.06212/applfast-cms-incjets-fnlo-arxiv-1512.06212.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-cms-incjets-fnlo-arxiv-1609.05331/applfast-cms-incjets-fnlo-arxiv-1609.05331.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-cms-incjets-fnlo-arxiv-2111.10431/applfast-cms-incjets-fnlo-arxiv-2111.10431.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-dijets-appl-arxiv-0010054/applfast-h1-dijets-appl-arxiv-0010054.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-dijets-appl-arxiv-0911.5678/applfast-h1-dijets-appl-arxiv-0911.5678.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-dijets-appl-arxiv-1406.4709/applfast-h1-dijets-appl-arxiv-1406.4709.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-dijets-appl-arxiv-1611.03421/applfast-h1-dijets-appl-arxiv-1611.03421.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-dijets-fnlo-arxiv-0010054/applfast-h1-dijets-fnlo-arxiv-0010054.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-dijets-fnlo-arxiv-0911.5678/applfast-h1-dijets-fnlo-arxiv-0911.5678.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-dijets-fnlo-arxiv-1406.4709/applfast-h1-dijets-fnlo-arxiv-1406.4709.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-dijets-fnlo-arxiv-1611.03421/applfast-h1-dijets-fnlo-arxiv-1611.03421.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-incjets-appl-arxiv-0010054/applfast-h1-incjets-appl-arxiv-0010054.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-incjets-appl-arxiv-0706.3722/applfast-h1-incjets-appl-arxiv-0706.3722.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-incjets-appl-arxiv-0911.5678/applfast-h1-incjets-appl-arxiv-0911.5678.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-incjets-appl-arxiv-1406.4709/applfast-h1-incjets-appl-arxiv-1406.4709.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-incjets-appl-arxiv-1611.03421/applfast-h1-incjets-appl-arxiv-1611.03421.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-incjets-fnlo-arxiv-0010054/applfast-h1-incjets-fnlo-arxiv-0010054.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-incjets-fnlo-arxiv-0706.3722/applfast-h1-incjets-fnlo-arxiv-0706.3722.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-incjets-fnlo-arxiv-0911.5678/applfast-h1-incjets-fnlo-arxiv-0911.5678.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-incjets-fnlo-arxiv-1406.4709/applfast-h1-incjets-fnlo-arxiv-1406.4709.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-h1-incjets-fnlo-arxiv-1611.03421/applfast-h1-incjets-fnlo-arxiv-1611.03421.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-zeus-dijets-appl-arxiv-1010.6167/applfast-zeus-dijets-appl-arxiv-1010.6167.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-zeus-dijets-fnlo-arxiv-1010.6167/applfast-zeus-dijets-fnlo-arxiv-1010.6167.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-zeus-incjets-appl-arxiv-0208037/applfast-zeus-incjets-appl-arxiv-0208037.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-zeus-incjets-appl-arxiv-0608048/applfast-zeus-incjets-appl-arxiv-0608048.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-zeus-incjets-fnlo-arxiv-0208037/applfast-zeus-incjets-fnlo-arxiv-0208037.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applfast/applfast-zeus-incjets-fnlo-arxiv-0608048/applfast-zeus-incjets-fnlo-arxiv-0608048.tgz"

    # Group: applgrid
    "https://ploughshare.web.cern.ch/ploughshare/db/applgrid/applgrid-cms-incjets-arxiv-1609.05331/applgrid-cms-incjets-arxiv-1609.05331.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/applgrid/applgrid-cms-top-173.3-arxiv-1703.01630/applgrid-cms-top-173.3-arxiv-1703.01630.tgz"

    # Group: atlas
    "https://ploughshare.web.cern.ch/ploughshare/db/atlas/atlas-atlas-dijets-13tev-arxiv-1711.02692/atlas-atlas-dijets-13tev-arxiv-1711.02692.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/atlas/atlas-atlas-dijets-arxiv-1112.6297/atlas-atlas-dijets-arxiv-1112.6297.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/atlas/atlas-atlas-dijets-arxiv-1312.3524/atlas-atlas-dijets-arxiv-1312.3524.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/atlas/atlas-atlas-incjet-13tev-arxiv-1711.02692/atlas-atlas-incjet-13tev-arxiv-1711.02692.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/atlas/atlas-atlas-incljets-8tev-arxiv-1706.03192/atlas-atlas-incljets-8tev-arxiv-1706.03192.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/atlas/atlas-atlas-incljets-arxiv-1009.5908v2/atlas-atlas-incljets-arxiv-1009.5908v2.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/atlas/atlas-atlas-incljets-arxiv-1112.6297/atlas-atlas-incljets-arxiv-1112.6297.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/atlas/atlas-atlas-incljets-arxiv-1304.4739/atlas-atlas-incljets-arxiv-1304.4739.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/atlas/atlas-atlas-incljets-arxiv-1410.8857/atlas-atlas-incljets-arxiv-1410.8857.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/atlas/atlas-atlas-wjets-arxiv-1711.03296/atlas-atlas-wjets-arxiv-1711.03296.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/atlas/atlas-atlas-wpm-arxiv-1109.5141/atlas-atlas-wpm-arxiv-1109.5141.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/atlas/atlas-atlas-wpm-arxiv-1612.03016/atlas-atlas-wpm-arxiv-1612.03016.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/atlas/atlas-atlas-z0-arxiv-1109.5141/atlas-atlas-z0-arxiv-1109.5141.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/atlas/atlas-atlas-z0-arxiv-1612.03016/atlas-atlas-z0-arxiv-1612.03016.tgz"

    # Group: fastnlo
    "https://ploughshare.web.cern.ch/ploughshare/db/fastnlo/fastnlo-cms-dijets-arxiv-1703.09986/fastnlo-cms-dijets-arxiv-1703.09986.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/fastnlo/fastnlo-cms-incjets-arxiv-1512.06212/fastnlo-cms-incjets-arxiv-1512.06212.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/fastnlo/fastnlo-cms-incjets-arxiv-1605.04436/fastnlo-cms-incjets-arxiv-1605.04436.tgz"

    # Group: xfitter
    "https://ploughshare.web.cern.ch/ploughshare/db/xfitter/xfitter-e615-dimuon-inspire-280845/xfitter-e615-dimuon-inspire-280845.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/xfitter/xfitter-na10-dimuon-194-inspire-212896/xfitter-na10-dimuon-194-inspire-212896.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/xfitter/xfitter-na10-dimuon-286-inspire-212896/xfitter-na10-dimuon-286-inspire-212896.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/xfitter/xfitter-wa70-pi+prompt-photon-inspire-250394/xfitter-wa70-pi+prompt-photon-inspire-250394.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/xfitter/xfitter-wa70-pi-prompt-photon-inspire-250394/xfitter-wa70-pi-prompt-photon-inspire-250394.tgz"

    ## Mitov's ttbar tables

    "http://www.precision.hep.phy.cam.ac.uk/wp-content/results/fastNLOtables/LHC8-ttbar-2dim-CMS.tar.gz"
    "http://www.precision.hep.phy.cam.ac.uk/wp-content/results/fastNLOtables/2d-ttbar-LHC-13TeV-3_masses.tar.gz"
    "http://www.precision.hep.phy.cam.ac.uk/wp-content/results/fastNLOtables/LHC13-ttbar-CMS_bin-fastNLO.tar.gz"
    "http://www.precision.hep.phy.cam.ac.uk/wp-content/results/fastNLOtables/LHC13-ttbar-CMS_bin-fastNLO-172_5.tar.gz"
    "http://www.precision.hep.phy.cam.ac.uk/wp-content/results/fastNLOtables/fastNLO-ttbar-NNLO-LHC8-173_3-bin1.tar.gz"
    "http://www.precision.hep.phy.cam.ac.uk/wp-content/results/fastNLOtables/fastNLO-ttbar-NNLO-LHC13-173_3-bin1.tar.gz"
)

tmp=$(mktemp -d)
cd "${tmp}"

trap "cd && rm -rf ${tmp}" EXIT SIGKILL

for grid in ${grids[@]}; do
    archive=${grid##*/}
    wget --no-verbose --no-clobber "${grid}" -O /tmp/"${archive}" || true
    mkdir subdir
    tar xzf /tmp/"${archive}" -C subdir

    for grid in $(find subdir \( -name '*.appl' -o -name '*.tab' -o -name '*.tab.gz' -o -name '*.root' \)); do
        if [[ $(basename $grid) =~ ^\..* ]]; then
            continue
        fi

        converted_grid="${grid}".pineappl.lz4
        pineappl import --accuracy 1e-12 "${grid}" "${converted_grid}" NNPDF31_nnlo_as_0118_luxqed
        du -h "${grid}" "${converted_grid}"
    done

    rm -r subdir
done
