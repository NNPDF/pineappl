#!/bin/bash -x

set -euo pipefail

grids=(
    ## Ploughshare

    # Group: pinejet
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-atlas-wm-arxiv-1109.5141/pinejet-atlas-wm-arxiv-1109.5141.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-atlas-wm-arxiv-1603.09222/pinejet-atlas-wm-arxiv-1603.09222.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-atlas-wm-arxiv-1612.03016/pinejet-atlas-wm-arxiv-1612.03016.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-atlas-wp-arxiv-1109.5141/pinejet-atlas-wp-arxiv-1109.5141.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-atlas-wp-arxiv-1603.09222/pinejet-atlas-wp-arxiv-1603.09222.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-atlas-wp-arxiv-1612.03016/pinejet-atlas-wp-arxiv-1612.03016.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-atlas-z0-arxiv-1109.5141/pinejet-atlas-z0-arxiv-1109.5141.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-atlas-z0-arxiv-1305.4192/pinejet-atlas-z0-arxiv-1305.4192.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-atlas-z0-arxiv-1404.1212/pinejet-atlas-z0-arxiv-1404.1212.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-atlas-z0-arxiv-1603.09222/pinejet-atlas-z0-arxiv-1603.09222.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-atlas-z0-arxiv-1606.01736/pinejet-atlas-z0-arxiv-1606.01736.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-atlas-z0-arxiv-1612.03016/pinejet-atlas-z0-arxiv-1612.03016.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-atlas-z0-arxiv-1710.05167/pinejet-atlas-z0-arxiv-1710.05167.tgz"
    #"https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-cdf-z0-arxiv-0908.3914/pinejet-cdf-z0-arxiv-0908.3914.tgz" # fails due to ppbar PDF mismatch
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-cms-wm-arxiv-1206.2598/pinejet-cms-wm-arxiv-1206.2598.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-cms-wm-arxiv-1312.6283/pinejet-cms-wm-arxiv-1312.6283.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-cms-wm-arxiv-1603.01803/pinejet-cms-wm-arxiv-1603.01803.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-cms-wp-arxiv-1206.2598/pinejet-cms-wp-arxiv-1206.2598.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-cms-wp-arxiv-1312.6283/pinejet-cms-wp-arxiv-1312.6283.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-cms-wp-arxiv-1603.01803/pinejet-cms-wp-arxiv-1603.01803.tgz"
    #"https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-cms-z0-arxiv-1310.7291/pinejet-cms-z0-arxiv-1310.7291.tgz" # fails due to static-scale optimization
    #"https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-d0-wm-arxiv-1309.2591/pinejet-d0-wm-arxiv-1309.2591.tgz" # fails du to ppbar PDF mismatch
    #"https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-d0-wp-arxiv-1309.2591/pinejet-d0-wp-arxiv-1309.2591.tgz" # fails du to ppbar PDF mismatch
    #"https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-d0-z0-arxiv-0702025/pinejet-d0-z0-arxiv-0702025.tgz" # fails du to ppbar PDF mismatch
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-lhcb-wm-arxiv-1505.07024/pinejet-lhcb-wm-arxiv-1505.07024.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-lhcb-wm-arxiv-1511.08039/pinejet-lhcb-wm-arxiv-1511.08039.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-lhcb-wp-arxiv-1505.07024/pinejet-lhcb-wp-arxiv-1505.07024.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-lhcb-wp-arxiv-1511.08039/pinejet-lhcb-wp-arxiv-1511.08039.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-lhcb-z0-arxiv-1505.07024/pinejet-lhcb-z0-arxiv-1505.07024.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-lhcb-z0-arxiv-1511.08039/pinejet-lhcb-z0-arxiv-1511.08039.tgz"
    "https://ploughshare.web.cern.ch/ploughshare/db/pinejet/pinejet-lhcb-z0-arxiv-1607.06495/pinejet-lhcb-z0-arxiv-1607.06495.tgz"
)

tmp=$(mktemp -d)
cd "${tmp}"

trap "cd && rm -rf ${tmp}" EXIT SIGKILL

for grid in ${grids[@]}; do
    archive=${grid##*/}
    wget --no-verbose --no-clobber "${grid}" -O /tmp/"${archive}" || true
    mkdir subdir
    tar xzf /tmp/"${archive}" -C subdir

    for grid in $(find subdir \( -name '*.lz4' \)); do
        if [[ $(basename $grid) =~ ^\..* ]]; then
            continue
        fi

        converted_grid="${grid}".appl
        reimported_grid="${grid}".new.pineappl.lz4
        pineappl export --accuracy 1e-12 "${grid}" "${converted_grid}" NNPDF31_nnlo_as_0118_luxqed
        pineappl import --accuracy 1e-12 "${converted_grid}" "${reimported_grid}" NNPDF31_nnlo_as_0118_luxqed
        du -h "${converted_grid}" "${reimported_grid}"
    done

    rm -r subdir
done
