#!/bin/bash

# exit script at the first sign of an error
set -o errexit

# the following exits if undeclared variables are used
set -o nounset

# exit if some program in a pipeline fails
set -o pipefail

# the directory that Madgraph5_aMC@NLO will use
dataset=NNPDF_DY_14TEV_40_PHENO

# make this number smaller to get more accurate results
accuracy=0.01

# the name of the grid
grid=${dataset}.pineappl

# try to find mg5_aMC
if ! which mg5_aMC 2> /dev/null; then
    echo "The binary \`mg5_aMC\` wasn't found. Add Madgraph5_aMC@NLO's bin folder to PATH." >&2
    exit 1
fi

# try to find pineappl
if ! which pineappl 2> /dev/null; then
    echo "The binary \`pineappl\` wasn't found. Please install it." >&2
    exit 1
fi

# step 1: generate code
cat > output.txt <<EOF
set complex_mass_scheme True
import model loop_qcd_qed_sm_Gmu
define p = p b b~
define j = p
generate p p > mu+ mu- [QCD QED]
output ${dataset}
quit
EOF

mg5_aMC output.txt

cd "${dataset}"

# step 2: implement custom scale choice
patch -p1 <<EOF
--- NLO/SubProcesses/setscales.f  2020-05-21 17:23:55.126143088 +0200
+++ NLO/SubProcesses/setscales.f.new  2020-05-21 17:27:26.262700419 +0200
@@ -527,6 +527,17 @@
       integer i,j
       character*80 temp_scale_id
       common/ctemp_scale_id/temp_scale_id
+      integer iPDG_reco(nexternal)
+      double precision ppl(0:3), pplb(0:3), ppv(0:3), xmll
+      double precision p_reco(0:4,nexternal), p_in(0:4,nexternal)
+c     les houches accord stuff to identify particles
+c
+      integer idup(nexternal,maxproc),mothup(2,nexternal,maxproc),
+     &        icolup(2,nexternal,maxflow),niprocs
+      common /c_leshouche_inc/idup,mothup,icolup,niprocs
+c Masses of external particles
+      double precision pmass(nexternal)
+      common/to_mass/pmass
 c
       tmp=0
       if(ickkw.eq.-1)then
@@ -568,10 +579,104 @@
 cc                 dynamical_scale_choice = 10                                   cc
 cc      in the run_card (run_card.dat)                                           cc
 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
-         write(*,*) "User-defined scale not set"
-         stop 1
-         temp_scale_id='User-defined dynamical scale' ! use a meaningful string
-         tmp = 0d0
+         temp_scale_id='Mll' ! use a meaningful string
+         tmp = -1d0
+         do i=1,nexternal
+           p_in(0:3,i) = pp(0:3,i)
+           p_in(4,i) = pmass(i)
+         enddo
+         call recombine_momenta(rphreco, etaphreco, lepphreco, quarkphreco,
+     $                          p_in, idup(1,1), p_reco, iPDG_reco)
+
+         do j = nincoming+1, nexternal
+           if (iPDG_reco(j).eq.13) ppl(0:3)=p_reco(0:3,j)
+           if (iPDG_reco(j).eq.-13) pplb(0:3)=p_reco(0:3,j)
+         enddo
+         do i=0,3
+           ppv(i)=ppl(i)+pplb(i)
+         enddo
+
+         xmll=sqrt(ppv(0)**2-ppv(1)**2-ppv(2)**2-ppv(3)**2)
+
+         if (xmll.lt.40d0) then
+           write (*,*) "error: event outside bins", xmll
+         else if (xmll.lt.45d0) then
+           tmp=0.5*(45d0+40d0)
+         else if (xmll.lt.50d0) then
+           tmp=0.5*(45d0+50d0)
+         else if (xmll.lt.55d0) then
+           tmp=0.5*(50d0+55d0)
+         else if (xmll.lt.60d0) then
+           tmp=0.5*(55d0+60d0)
+         else if (xmll.lt.64d0) then
+           tmp=0.5*(60d0+64d0)
+         else if (xmll.lt.68d0) then
+           tmp=0.5*(64d0+68d0)
+         else if (xmll.lt.72d0) then
+           tmp=0.5*(68d0+72d0)
+         else if (xmll.lt.76d0) then
+           tmp=0.5*(72d0+76d0)
+         else if (xmll.lt.81d0) then
+           tmp=0.5*(76d0+81d0)
+         else if (xmll.lt.86d0) then
+           tmp=0.5*(81d0+86d0)
+         else if (xmll.lt.91d0) then
+           tmp=0.5*(86d0+91d0)
+         else if (xmll.lt.96d0) then
+           tmp=0.5*(91d0+96d0)
+         else if (xmll.lt.101d0) then
+           tmp=0.5*(96d0+101d0)
+         else if (xmll.lt.106d0) then
+           tmp=0.5*(101d0+106d0)
+         else if (xmll.lt.110d0) then
+           tmp=0.5*(106d0+110d0)
+         else if (xmll.lt.115d0) then
+           tmp=0.5*(110d0+115d0)
+         else if (xmll.lt.120d0) then
+           tmp=0.5*(115d0+120d0)
+         else if (xmll.lt.126d0) then
+           tmp=0.5*(120d0+126d0)
+         else if (xmll.lt.133d0) then
+           tmp=0.5*(126d0+133d0)
+         else if (xmll.lt.141d0) then
+           tmp=0.5*(133d0+141d0)
+         else if (xmll.lt.150d0) then
+           tmp=0.5*(141d0+150d0)
+         else if (xmll.lt.160d0) then
+           tmp=0.5*(150d0+160d0)
+         else if (xmll.lt.171d0) then
+           tmp=0.5*(160d0+171d0)
+         else if (xmll.lt.185d0) then
+           tmp=0.5*(171d0+185d0)
+         else if (xmll.lt.200d0) then
+           tmp=0.5*(185d0+200d0)
+         else if (xmll.lt.220d0) then
+           tmp=0.5*(200d0+220d0)
+         else if (xmll.lt.243d0) then
+           tmp=0.5*(220d0+243d0)
+         else if (xmll.lt.273d0) then
+           tmp=0.5*(243d0+273d0)
+         else if (xmll.lt.320d0) then
+           tmp=0.5*(273d0+320d0)
+         else if (xmll.lt.380d0) then
+           tmp=0.5*(320d0+380d0)
+         else if (xmll.lt.440d0) then
+           tmp=0.5*(380d0+440d0)
+         else if (xmll.lt.510d0) then
+           tmp=0.5*(440d0+510d0)
+         else if (xmll.lt.600d0) then
+           tmp=0.5*(510d0+600d0)
+         else if (xmll.lt.700d0) then
+           tmp=0.5*(600d0+700d0)
+         else if (xmll.lt.830d0) then
+           tmp=0.5*(700d0+830d0)
+         else if (xmll.lt.1000d0) then
+           tmp=0.5*(830d0+1000d0)
+         else if (xmll.lt.1500d0) then
+           tmp=0.5*(1000d0+1500d0)
+         else if (xmll.lt.3000d0) then
+           tmp=0.5*(1500d0+3000d0)
+         endif
 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 cc      USER-DEFINED SCALE: END OF USER CODE                                     cc
 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
EOF

# step 3: write analysis file
cat > FixedOrderAnalysis/"${dataset}".f <<EOF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_begin(nwgt,weights_info)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      integer nwgt
      character*(*) weights_info(*)

      call set_error_estimation(1)
      call HwU_inithist(nwgt,weights_info)
      call HwU_book(1,'mll', 38, 0d0, 38d0)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_end(dummy)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      double precision dummy
      call HwU_write_file
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_fill(p,istatus,ipdg,wgts,ibody)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'nexternal.inc'
      include 'cuts.inc'
      integer istatus(nexternal)
      integer iPDG(nexternal)
      integer ibody
      double precision p(0:4,nexternal)
      double precision wgts(*)
      double precision ppl(0:3), pplb(0:3), ppv(0:3), xmll
      double precision obs
      integer i

      double precision p_reco(0:4,nexternal)
      integer iPDG_reco(nexternal),grid

      call recombine_momenta(rphreco, etaphreco, lepphreco, quarkphreco,
     $                       p, iPDG, p_reco, iPDG_reco)

      do i = nincoming+1, nexternal
        if (iPDG_reco(i).eq.13) ppl(0:3)=p_reco(0:3,i)
        if (iPDG_reco(i).eq.-13) pplb(0:3)=p_reco(0:3,i)
      enddo
      do i=0,3
        ppv(i)=ppl(i)+pplb(i)
      enddo

      xmll=sqrt(ppv(0)**2-ppv(1)**2-ppv(2)**2-ppv(3)**2)
      obs=-1d0

      if (xmll.lt.40d0) then
        obs=-2
      else if (xmll.lt.45d0) then
        obs=0.5d0
      else if (xmll.lt.50d0) then
        obs=1.5d0
      else if (xmll.lt.55d0) then
        obs=2.5d0
      else if (xmll.lt.60d0) then
        obs=3.5d0
      else if (xmll.lt.64d0) then
        obs=4.5d0
      else if (xmll.lt.68d0) then
        obs=5.5d0
      else if (xmll.lt.72d0) then
        obs=6.5d0
      else if (xmll.lt.76d0) then
        obs=7.5d0
      else if (xmll.lt.81d0) then
        obs=8.5d0
      else if (xmll.lt.86d0) then
        obs=9.5d0
      else if (xmll.lt.91d0) then
        obs=10.5d0
      else if (xmll.lt.96d0) then
        obs=11.5d0
      else if (xmll.lt.101d0) then
        obs=12.5d0
      else if (xmll.lt.106d0) then
        obs=13.5d0
      else if (xmll.lt.110d0) then
        obs=14.5d0
      else if (xmll.lt.115d0) then
        obs=15.5d0
      else if (xmll.lt.120d0) then
        obs=16.5d0
      else if (xmll.lt.126d0) then
        obs=17.5d0
      else if (xmll.lt.133d0) then
        obs=18.5d0
      else if (xmll.lt.141d0) then
        obs=19.5d0
      else if (xmll.lt.150d0) then
        obs=20.5d0
      else if (xmll.lt.160d0) then
        obs=21.5d0
      else if (xmll.lt.171d0) then
        obs=22.5d0
      else if (xmll.lt.185d0) then
        obs=23.5d0
      else if (xmll.lt.200d0) then
        obs=24.5d0
      else if (xmll.lt.220d0) then
        obs=25.5d0
      else if (xmll.lt.243d0) then
        obs=26.5d0
      else if (xmll.lt.273d0) then
        obs=27.5d0
      else if (xmll.lt.320d0) then
        obs=28.5d0
      else if (xmll.lt.380d0) then
        obs=29.5d0
      else if (xmll.lt.440d0) then
        obs=30.5d0
      else if (xmll.lt.510d0) then
        obs=31.5d0
      else if (xmll.lt.600d0) then
        obs=32.5d0
      else if (xmll.lt.700d0) then
        obs=33.5d0
      else if (xmll.lt.830d0) then
        obs=34.5d0
      else if (xmll.lt.1000d0) then
        obs=35.5d0
      else if (xmll.lt.1500d0) then
        obs=36.5d0
      else if (xmll.lt.3000d0) then
        obs=37.5d0
      else
        obs=-3
      endif

      if (obs.gt.0d0) then
        call HwU_fill(1,obs,wgts)
      endif

 999  return
      end
EOF

cd -

# step 4: activate analysis file
sed -i "s/analysis_HwU_template/${dataset}/g" "${dataset}"/Cards/FO_analyse_card.dat

# step 5: generate the PineAPPL file
cat > launch.txt <<EOF
launch NNPDF_DY_14TEV_40_PHENO
fixed_order = ON
set maxjetflavor 5
set gf 1.1663787e-5
set mh 125.0
set mt 172.5
set mw 80.352
set mz 91.1535
set wh 4.07468e-3
set wt 1.37758
set ww 2.084
set wz 2.4943
set ebeam1 7000
set ebeam2 7000
set pdlabel lhapdf
set lhaid 324900
set dynamical_scale_choice 10
set reweight_scale True
set ptl = 15.0
set etal = 2.4
set mll_sf = 40.0
set rphreco = 0.1
set req_acc_FO ${accuracy}
# the following is the only switch needed to activate PineAPPL
set pineappl True
done
quit
EOF

mg5_aMC launch.txt

# step 6: copy the PineAPPL file
cp ${dataset}/Events/run_01/amcblast_obs_0.pineappl "${grid}"

# step 7: change the bin limits from 0,1,2,3,... to the actual values
pineappl remap "${grid}" "${grid}".tmp '40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000,1500,3000'
mv "${grid}".tmp "${grid}"

# step 8: add some metadata (used mainly by the plot script)
pineappl set "${grid}" "${grid}".tmp \
    --entry description 'Differential Drell–Yan cross section at 14 TeV' \
    --entry x1_label 'Mll' \
    --entry x1_label_tex '$M_{\ell\bar{\ell}}$' \
    --entry x1_unit 'GeV' \
    --entry y_label 'dsig/dMll' \
    --entry y_label_tex '$\frac{\mathrm{d}\sigma}{\mathrm{d}M_{\ell\bar{\ell}}}$' \
    --entry y_unit 'pb/GeV'
mv "${grid}".tmp "${grid}"

# step 9: compress the grid if we find `lz4`
lz4=$(which lz4 2> /dev/null || true)

if [[ -x ${lz4} ]]; then
    lz4 -9 "${grid}"
    rm "${grid}"
    grid="${grid}".lz4
fi

cat <<EOF
Generated ${grid}.

Try using:
  - pineappl convolute ${grid} LHAPDF_SET_NAME
  - pineappl --silence-lhapdf plot ${grid} LHAPDF_SET_NAME1 LHAPDF_SET_NAME2 ... > plot_script.py
  - pineappl --help
EOF
