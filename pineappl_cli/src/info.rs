use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{ArgGroup, Parser, ValueHint};
use itertools::Itertools;
use std::path::PathBuf;

/// Shows information about the grid.
#[derive(Parser)]
#[clap(group = ArgGroup::new("mode").required(true), name = "info")]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// For each order print a list of the largest EW order.
    #[clap(group = "mode", long)]
    ew: bool,
    /// Gets an internal key-value pair.
    #[clap(group = "mode", long, takes_value = true, value_name = "key")]
    get: Option<String>,
    /// Show all keys stored in the grid.
    #[clap(group = "mode", long)]
    keys: bool,
    /// For each order print a list of the largest QCD order.
    #[clap(group = "mode", long)]
    qcd: bool,
    /// Shows all key-value pairs stored in the grid.
    #[clap(group = "mode", long)]
    show: bool,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let mut grid = helpers::read_grid(&self.input)?;

        if self.ew || self.qcd {
            let mut sorted_grid_orders: Vec<_> = grid
                .orders()
                .iter()
                .filter(|order| (order.logxir == 0) && (order.logxif == 0))
                .collect();
            sorted_grid_orders.sort();

            let orders = sorted_grid_orders
                .into_iter()
                .group_by(|order| order.alphas + order.alpha)
                .into_iter()
                .map(|mut iter| {
                    if self.qcd {
                        iter.1.next().unwrap()
                    } else {
                        iter.1.last().unwrap()
                    }
                })
                .map(|order| {
                    if order.alphas == 0 {
                        format!("a{}", order.alpha)
                    } else if order.alpha == 0 {
                        format!("as{}", order.alphas)
                    } else {
                        format!("as{}a{}", order.alphas, order.alpha)
                    }
                })
                .collect::<Vec<_>>()
                .join(",");

            println!("{}", orders);
        } else if let Some(ref key) = self.get {
            grid.upgrade();

            if let Some(key_values) = grid.key_values() {
                if let Some(value) = key_values.get(key) {
                    println!("{}", value);
                }
            } else {
                unreachable!();
            }
        } else if self.keys {
            grid.upgrade();

            if let Some(key_values) = grid.key_values() {
                let mut vector = key_values.iter().collect::<Vec<_>>();
                vector.sort();

                for (key, _) in &vector {
                    println!("{}", key);
                }
            } else {
                unreachable!();
            }
        } else if self.show {
            grid.upgrade();

            if let Some(key_values) = grid.key_values() {
                let mut vector = key_values.iter().collect::<Vec<_>>();
                vector.sort();

                for (key, value) in &vector {
                    println!("{}: {}", key, value);
                }
            } else {
                unreachable!();
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl-info 

Shows information about the grid

USAGE:
    pineappl info <--ew|--get <key>|--keys|--qcd|--show> <INPUT>

ARGS:
    <INPUT>    Path to the input grid

OPTIONS:
        --ew           For each order print a list of the largest EW order
        --get <key>    Gets an internal key-value pair
    -h, --help         Print help information
        --keys         Show all keys stored in the grid
        --qcd          For each order print a list of the largest QCD order
        --show         Shows all key-value pairs stored in the grid
";

    const KEYS_STR: &str = "arxiv
description
hepdata
initial_state_1
initial_state_2
lumi_id_types
mg5amc_repo
mg5amc_revno
nnpdf_id
pineappl_gitversion
results
runcard
runcard_gitversion
x1_label
x1_label_tex
x1_unit
y_label
y_label_tex
y_unit
";

    const SHOW_STR: &str = r#"arxiv: 1505.07024
description: LHCb differential W-boson production cross section at 7 TeV
hepdata: 10.17182/hepdata.2114.v1/t4
initial_state_1: 2212
initial_state_2: 2212
lumi_id_types: pdg_mc_ids
mg5amc_repo: 
mg5amc_revno: 
nnpdf_id: LHCBWZMU7TEV
pineappl_gitversion: v0.4.1-36-gdbdb5d0
results: ----------------------------------------------------------------------
   PineAPPL         MC        sigma      central         min      max 
                              1/100   sigma   1/1000   1/1000   1/1000
----------------------------------------------------------------------
 1.876381e+02  1.876313e+02   0.082   0.044   0.0360   0.0396   0.0317
 1.726078e+02  1.726041e+02   0.082   0.026   0.0211   0.0246   0.0166
 1.500070e+02  1.500056e+02   0.051   0.019   0.0095   0.0103   0.0075
 1.212883e+02  1.212890e+02   0.047   0.012   0.0056   0.0052   0.0068
 9.046672e+01  9.046795e+01   0.057   0.024   0.0136   0.0127   0.0146
 6.145558e+01  6.145650e+01   0.063   0.024   0.0151   0.0116   0.0193
 5.785102e+01  5.784687e+01   0.075   0.095   0.0717   0.0874   0.0518
 1.377203e+01  1.376219e+01   0.119   0.599   0.7153   0.7646   0.6465

runcard: <LesHouchesEvents version="3.0">
<header>
<!--
#*********************************************************************
#                                                                    *
#                        MadGraph5_aMC@NLO                           *
#                                                                    *
#                           Going Beyond                             *
#                                                                    *
#                   http://madgraph.hep.uiuc.edu                     *
#                   http://madgraph.phys.ucl.ac.be                   *
#                   http://amcatnlo.cern.ch                          *
#                                                                    *
#                     The MadGraph5_aMC@NLO team                     *
#                                                                    *
#....................................................................*
#                                                                    *
# This file contains all the information necessary to reproduce      *
# the events generated:                                              *
#                                                                    *
# 1. software version                                                *
# 2. proc_card          : code generation info including model       *
# 3. param_card         : model primary parameters in the LH format  *
# 4. run_card           : running parameters (collider and cuts)     *
# 5. pythia_card        : present only if pythia has been run        *
# 6. pgs_card           : present only if pgs has been run           *
# 7. delphes_cards      : present only if delphes has been run       *
#                                                                    *
#                                                                    *
#*********************************************************************
-->
<MGVersion>
3.1.0
</MGVersion>
<MGRunCard>
<![CDATA[
#***********************************************************************
#                        MadGraph5_aMC@NLO                             *
#                                                                      *
#                      run_card.dat aMC@NLO                            *
#                                                                      *
#  This file is used to set the parameters of the run.                 *
#                                                                      *
#  Some notation/conventions:                                          *
#                                                                      *
#   Lines starting with a hash (#) are info or comments                *
#                                                                      *
#   mind the format:   value    = variable     ! comment               *
#                                                                      *
#   Some of the values of variables can be list. These can either be   *
#   comma or space separated.                                          *
#                                                                      *
#   To display additional parameter, you can use the command:          *
#      update to_full                                                  *
#***********************************************************************
#
#*******************                                                 
# Running parameters
#*******************                                                 
#
#***********************************************************************
# Tag name for the run (one word)                                      *
#***********************************************************************
  tag_1	= run_tag ! name of the run 
#***********************************************************************
# Number of LHE events (and their normalization) and the required      *
# (relative) accuracy on the Xsec.                                     *
# These values are ignored for fixed order runs                        *
#***********************************************************************
  10000	= nevents ! Number of unweighted events requested 
  -1.0	= req_acc ! Required accuracy (-1=auto determined from nevents)
  -1	= nevt_job ! Max number of events per job in event generation. 
                 !  (-1= no split).
#***********************************************************************
# Normalize the weights of LHE events such that they sum or average to *
# the total cross section                                              *
#***********************************************************************
  average	= event_norm ! valid settings: average, sum, bias
#***********************************************************************
# Number of points per itegration channel (ignored for aMC@NLO runs)   *
#***********************************************************************
  0.0002	= req_acc_fo ! Required accuracy (-1=ignored, and use the 
 	                   ! number of points and iter. below)
# These numbers are ignored except if req_acc_FO is equal to -1
  5000	= npoints_fo_grid ! number of points to setup grids
  4	= niters_fo_grid ! number of iter. to setup grids
  10000	= npoints_fo ! number of points to compute Xsec
  6	= niters_fo ! number of iter. to compute Xsec
#***********************************************************************
# Random number seed                                                   *
#***********************************************************************
  0	= iseed ! rnd seed (0=assigned automatically=default))
#***********************************************************************
# Collider type and energy                                             *
#***********************************************************************
  1	= lpp1 ! beam 1 type (0 = no PDF)
  1	= lpp2 ! beam 2 type (0 = no PDF)
  3500.0	= ebeam1 ! beam 1 energy in GeV
  3500.0	= ebeam2 ! beam 2 energy in GeV
#***********************************************************************
# PDF choice: this automatically fixes also alpha_s(MZ) and its evol.  *
#***********************************************************************
  lhapdf	= pdlabel ! PDF set
  324900	= lhaid ! If pdlabel=lhapdf, this is the lhapdf number. Only 
              ! numbers for central PDF sets are allowed. Can be a list; 
              ! PDF sets beyond the first are included via reweighting.
#***********************************************************************
# Include the NLO Monte Carlo subtr. terms for the following parton    *
# shower (HERWIG6 | HERWIGPP | PYTHIA6Q | PYTHIA6PT | PYTHIA8)         *
# WARNING: PYTHIA6PT works only for processes without FSR!!!!          *
#***********************************************************************
  HERWIG6	= parton_shower 
  1.0	= shower_scale_factor ! multiply default shower starting
                                  ! scale by this factor
#***********************************************************************
# Renormalization and factorization scales                             *
# (Default functional form for the non-fixed scales is the sum of      *
# the transverse masses divided by two of all final state particles    * 
# and partons. This can be changed in SubProcesses/set_scales.f or via *
# dynamical_scale_choice option)                                       *
#***********************************************************************
  True	= fixed_ren_scale ! if .true. use fixed ren scale
  True	= fixed_fac_scale ! if .true. use fixed fac scale
  80.352	= mur_ref_fixed ! fixed ren reference scale 
  80.352	= muf_ref_fixed ! fixed fact reference scale
  -1	= dynamical_scale_choice ! Choose one (or more) of the predefined
           ! dynamical choices. Can be a list; scale choices beyond the
           ! first are included via reweighting
  1.0	= mur_over_ref ! ratio of current muR over reference muR
  1.0	= muf_over_ref ! ratio of current muF over reference muF
#*********************************************************************** 
# Reweight variables for scale dependence and PDF uncertainty          *
#***********************************************************************
  1.0, 2.0, 0.5	= rw_rscale ! muR factors to be included by reweighting
  1.0, 2.0, 0.5	= rw_fscale ! muF factors to be included by reweighting
  True	= reweight_scale ! Reweight to get scale variation using the 
            ! rw_rscale and rw_fscale factors. Should be a list of 
            ! booleans of equal length to dynamical_scale_choice to
            ! specify for which choice to include scale dependence.
  False	= reweight_pdf ! Reweight to get PDF uncertainty. Should be a
            ! list booleans of equal length to lhaid to specify for
            !  which PDF set to include the uncertainties.
#***********************************************************************
# Store reweight information in the LHE file for off-line model-       *
# parameter reweighting at NLO+PS accuracy                             *
#***********************************************************************
  False	= store_rwgt_info ! Store info for reweighting in LHE file
#***********************************************************************
# ickkw parameter:                                                     *
#   0: No merging                                                      *
#   3: FxFx Merging - WARNING! Applies merging only at the hard-event  *
#      level. After showering an MLM-type merging should be applied as *
#      well. See http://amcatnlo.cern.ch/FxFx_merging.htm for details. *
#   4: UNLOPS merging (with pythia8 only). No interface from within    *
#      MG5_aMC available, but available in Pythia8.                    *
#  -1: NNLL+NLO jet-veto computation. See arxiv:1412.8408 [hep-ph].    *
#***********************************************************************
  0	= ickkw 
#***********************************************************************
#
#***********************************************************************
# BW cutoff (M+/-bwcutoff*Gamma). Determines which resonances are      *
# written in the LHE event file                                        *
#***********************************************************************
  15.0	= bwcutoff 
#***********************************************************************
# Cuts on the jets. Jet clustering is performed by FastJet.            *
#  - If gamma_is_j, photons are also clustered                            *
#  - When matching to a parton shower, these generation cuts should be *
#    considerably softer than the analysis cuts.                       *
#  - More specific cuts can be specified in SubProcesses/cuts.f        *
#***********************************************************************
  1.0	= jetalgo ! FastJet jet algorithm (1=kT, 0=C/A, -1=anti-kT)
  0.7	= jetradius ! The radius parameter for the jet algorithm
  10.0	= ptj ! Min jet transverse momentum
  -1.0	= etaj ! Max jet abs(pseudo-rap) (a value .lt.0 means no cut)
  True	= gamma_is_j ! Wether to cluster photons as jets or not
#***********************************************************************
# Cuts on the charged leptons (e+, e-, mu+, mu-, tau+ and tau-)        *
# More specific cuts can be specified in SubProcesses/cuts.f           *
#***********************************************************************
  20.0	= ptl ! Min lepton transverse momentum
  4.5	= etal ! Max lepton abs(pseudo-rap) (a value .lt.0 means no cut)
  0.0	= drll ! Min distance between opposite sign lepton pairs
  0.0	= drll_sf ! Min distance between opp. sign same-flavor lepton pairs
  0.0	= mll ! Min inv. mass of all opposite sign lepton pairs
  30.0	= mll_sf ! Min inv. mass of all opp. sign same-flavor lepton pairs
#***********************************************************************
# Fermion-photon recombination parameters                              *
# If Rphreco=0, no recombination is performed                          *
#***********************************************************************
  0.1	= rphreco ! Minimum fermion-photon distance for recombination
  -1.0	= etaphreco ! Maximum abs(pseudo-rap) for photons to be recombined (a value .lt.0 means no cut)
  True	= lepphreco ! Recombine photons and leptons together
  True	= quarkphreco ! Recombine photons and quarks together
#***********************************************************************
# Photon-isolation cuts, according to hep-ph/9801442                   *
# Not applied if gamma_is_j                                            *
# When ptgmin=0, all the other parameters are ignored                  *
# More specific cuts can be specified in SubProcesses/cuts.f           *
#***********************************************************************
  20.0	= ptgmin ! Min photon transverse momentum
  -1.0	= etagamma ! Max photon abs(pseudo-rap)
  0.4	= r0gamma ! Radius of isolation code
  1.0	= xn ! n parameter of eq.(3.4) in hep-ph/9801442
  1.0	= epsgamma ! epsilon_gamma parameter of eq.(3.4) in hep-ph/9801442
  True	= isoem ! isolate photons from EM energy (photons and leptons)
#***********************************************************************
# Cuts associated to MASSIVE particles identified by their PDG codes.  *
# All cuts are applied to both particles and anti-particles, so use    *
# POSITIVE PDG CODES only. Example of the syntax is {6 : 100} or       *
# {6:100, 25:200} for multiple particles                               *
#***********************************************************************
  {}	= pt_min_pdg ! Min pT for a massive particle
  {}	= pt_max_pdg ! Max pT for a massive particle
  {}	= mxx_min_pdg ! inv. mass for any pair of (anti)particles
#***********************************************************************
# Use PineAPPL to generate PDF-independent fast-interpolation grid     *
# (https://zenodo.org/record/3992765#.X2EWy5MzbVo)                     *
#***********************************************************************
  True	= pineappl ! PineAPPL switch 
#***********************************************************************
#********************************************************************* 
#  Additional hidden parameters
#*********************************************************************
  5	= maxjetflavor # hidden_parameter
]]>
</MGRunCard>
<slha>
######################################################################
## PARAM_CARD AUTOMATICALY GENERATED BY MG5                       ####
######################################################################
###################################
## INFORMATION FOR MASS
###################################
BLOCK MASS # 
      6 1.725000e+02 # mt
      23 9.115350e+01 # mz
      24 8.035200e+01 # mw
      25 1.250000e+02 # mh
      1 0.000000e+00 # d : 0.0
      2 0.000000e+00 # u : 0.0
      3 0.000000e+00 # s : 0.0
      4 0.000000e+00 # c : 0.0
      5 0.000000e+00 # b : 0.0
      11 0.000000e+00 # e- : 0.0
      12 0.000000e+00 # ve : 0.0
      13 0.000000e+00 # mu- : 0.0
      14 0.000000e+00 # vm : 0.0
      15 0.000000e+00 # ta- : 0.0
      16 0.000000e+00 # vt : 0.0
      21 0.000000e+00 # g : 0.0
      22 0.000000e+00 # a : 0.0
      9000002 9.118800e+01 # ghz : mz
      9000003 8.041900e+01 # ghwp : mw
      9000004 8.041900e+01 # ghwm : mw
      250 9.118800e+01 # g0 : mz
      251 8.041900e+01 # g+ : mw
###################################
## INFORMATION FOR SMINPUTS
###################################
BLOCK SMINPUTS # 
      2 1.166379e-05 # gf
      3 1.180000e-01 # as (note that parameter not used if you use a pdf set)
###################################
## INFORMATION FOR DECAY
###################################
DECAY 6 1.377580e+00 # wt
DECAY 23 2.494300e+00 # wz
DECAY 24 2.084000e+00 # ww
DECAY 25 4.074680e-03 # wh
DECAY 1 0.000000e+00 # d : 0.0
DECAY 2 0.000000e+00 # u : 0.0
DECAY 3 0.000000e+00 # s : 0.0
DECAY 4 0.000000e+00 # c : 0.0
DECAY 5 0.000000e+00 # b : 0.0
DECAY 11 0.000000e+00 # e- : 0.0
DECAY 12 0.000000e+00 # ve : 0.0
DECAY 13 0.000000e+00 # mu- : 0.0
DECAY 14 0.000000e+00 # vm : 0.0
DECAY 15 0.000000e+00 # ta- : 0.0
DECAY 16 0.000000e+00 # vt : 0.0
DECAY 21 0.000000e+00 # g : 0.0
DECAY 22 0.000000e+00 # a : 0.0
DECAY 250 2.504790e+00 # g0 : wz
DECAY 251 2.092910e+00 # g+ : ww
###################################
## INFORMATION FOR QNUMBERS 9000001
###################################
BLOCK QNUMBERS 9000001 #  gha
      1 0 # 3 times electric charge
      2 1 # number of spin states (2s+1)
      3 1 # colour rep (1: singlet, 3: triplet, 8: octet)
      4 1 # particle/antiparticle distinction (0=own anti)
###################################
## INFORMATION FOR QNUMBERS 9000002
###################################
BLOCK QNUMBERS 9000002 #  ghz
      1 0 # 3 times electric charge
      2 1 # number of spin states (2s+1)
      3 1 # colour rep (1: singlet, 3: triplet, 8: octet)
      4 1 # particle/antiparticle distinction (0=own anti)
###################################
## INFORMATION FOR QNUMBERS 9000003
###################################
BLOCK QNUMBERS 9000003 #  ghwp
      1 3 # 3 times electric charge
      2 1 # number of spin states (2s+1)
      3 1 # colour rep (1: singlet, 3: triplet, 8: octet)
      4 1 # particle/antiparticle distinction (0=own anti)
###################################
## INFORMATION FOR QNUMBERS 9000004
###################################
BLOCK QNUMBERS 9000004 #  ghwm
      1 -3 # 3 times electric charge
      2 1 # number of spin states (2s+1)
      3 1 # colour rep (1: singlet, 3: triplet, 8: octet)
      4 1 # particle/antiparticle distinction (0=own anti)
###################################
## INFORMATION FOR QNUMBERS 9000005
###################################
BLOCK QNUMBERS 9000005 #  ghg
      1 0 # 3 times electric charge
      2 1 # number of spin states (2s+1)
      3 8 # colour rep (1: singlet, 3: triplet, 8: octet)
      4 1 # particle/antiparticle distinction (0=own anti)
###################################
## INFORMATION FOR QNUMBERS 250
###################################
BLOCK QNUMBERS 250 #  g0
      1 0 # 3 times electric charge
      2 1 # number of spin states (2s+1)
      3 1 # colour rep (1: singlet, 3: triplet, 8: octet)
      4 0 # particle/antiparticle distinction (0=own anti)
###################################
## INFORMATION FOR QNUMBERS 251
###################################
BLOCK QNUMBERS 251 #  g+
      1 3 # 3 times electric charge
      2 1 # number of spin states (2s+1)
      3 1 # colour rep (1: singlet, 3: triplet, 8: octet)
      4 1 # particle/antiparticle distinction (0=own anti)
</slha>
<run_settings>
order = NLO
fixed_order = ON
shower = OFF
madspin = OFF
reweight = OFF
madanalysis = OFF
runshower = False
</run_settings>
<foanalyse>
<![CDATA[
#######################################################################
#                             
# This file contains the settings for analyses to be linked to fixed
# order runs. Analysis files are meant to be put (or linked) inside
# <PROCDIR>/FixedOrderAnalysis/ (<PROCDIR> is the name of the exported
# process directory). See the
# <PROCDIR>/FixedOrderAnalysis/analysis_*_template.f file for details
# on how to write your own analysis.
#                                                                               
#######################################################################
#
# Analysis format.
# Can either be 'topdrawer', 'root', 'HwU', 'LHE' or 'none'.
# When choosing HwU, it comes with a GnuPlot wrapper. When choosing
# topdrawer, the histogramming package 'dbook.f' is included in the
# code, while when choosing root the 'rbook_fe8.f' and 'rbook_be8.cc'
# are included. If 'none' is chosen, all the other entries below have
# to be set empty.
FO_ANALYSIS_FORMAT = HwU
#
#
# Needed extra-libraries (FastJet is already linked):
FO_EXTRALIBS = 
#
# (Absolute) path to the extra libraries. Directory names should be
# separated by white spaces.
FO_EXTRAPATHS =
#
# (Absolute) path to the dirs containing header files needed by the
# libraries (e.g. C++ header files):
FO_INCLUDEPATHS =                      
#
# User's analysis (to be put in the <PROCDIR>/FixedOrderAnalysis/
# directory). Please use .o as extension and white spaces to separate
# files.
FO_ANALYSE = LHCB_WP_7TEV.o
#
#
## When linking with root, the following settings are a working
## example on lxplus (CERN) as of July 2014. When using this, comment
## out the lines above and replace <PATH_TO_ROOT> with the physical
## path to root,
## e.g. /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.11/x86_64-slc6-gcc46-dbg/root/
#FO_ANALYSIS_FORMAT = root
#FO_EXTRALIBS = Core Cint Hist Matrix MathCore RIO dl Thread
#FO_EXTRAPATHS = <PATH_TO_ROOT>/lib
#FO_INCLUDEPATHS = <PATH_TO_ROOT>/include
#FO_ANALYSE = analysis_root_template.o
]]>
</foanalyse>
</header>
</LesHouchesEvents>

runcard_gitversion: 7b42083
x1_label: etal
x1_label_tex: $\eta_{\bar{\ell}}$
x1_unit: 
y_label: disg/detal
y_label_tex: $\frac{\mathrm{d}\sigma}{\mathrm{d}\eta_{\bar{\ell}}}$
y_unit: pb
"#;

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["info", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }

    #[test]
    fn ew() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["info", "--ew", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout("a2,a3\n");
    }

    #[test]
    fn get_arxiv() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["info", "--get=arxiv", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout("1505.07024\n");
    }

    #[test]
    fn keys() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["info", "--keys", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout(KEYS_STR);
    }

    #[test]
    fn qcd() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["info", "--qcd", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout("a2,as1a2\n");
    }

    #[test]
    fn show() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["info", "--show", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout(SHOW_STR);
    }
}
