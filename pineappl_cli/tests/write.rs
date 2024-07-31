use assert_cmd::Command;
use assert_fs::{fixture::FileWriteStr, NamedTempFile};

const HELP_STR: &str = "Write a grid modified by various operations

Usage: pineappl write [OPTIONS] <INPUT> <OUTPUT>

Arguments:
  <INPUT>   Path to the input grid
  <OUTPUT>  Path of the modified PineAPPL file

Options:
      --cc1[=<ENABLE>]                 Charge conjugate the first initial state [possible values: true, false]
      --cc2[=<ENABLE>]                 Charge conjugate the second initial state [possible values: true, false]
      --dedup-channels[=<ULPS>]        Deduplicate channels assuming numbers differing by ULPS are the same
      --delete-bins <BIN1-BIN2,...>    Delete bins with the specified indices
      --delete-channels <CH1-CH2,...>  Delete channels with the specified indices
      --delete-key <KEY>               Delete an internal key-value pair
      --merge-bins <BIN1-BIN2,...>     Merge specific bins together
      --optimize[=<ENABLE>]            Optimize internal data structure to minimize memory and disk usage [possible values: true, false]
      --optimize-fk-table <OPTIMI>     Optimize internal data structure of an FkTable to minimize memory and disk usage [possible values: Nf6Ind, Nf6Sym, Nf5Ind, Nf5Sym, Nf4Ind, Nf4Sym, Nf3Ind, Nf3Sym]
      --remap <REMAPPING>              Modify the bin dimensions and widths
      --remap-norm <NORM>              Modify the bin normalizations with a common factor
      --remap-norm-ignore <DIM1,...>   Modify the bin normalizations by multiplying with the bin lengths for the given dimensions
      --rewrite-channel <IDX> <CHAN>   Rewrite the definition of the channel with index IDX
      --rewrite-order <IDX> <ORDER>    Rewrite the definition of the order with index IDX
      --rotate-pid-basis <BASIS>       Rotate the PID basis for this grid [possible values: PDG, EVOL]
  -s, --scale <SCALE>                  Scales all grids with the given factor
      --scale-by-bin <BIN1,BIN2,...>   Scale each bin with a different factor
      --scale-by-order <AS,AL,LR,LF>   Scales all grids with order-dependent factors
      --set-key-value <KEY> <VALUE>    Set an internal key-value pair
      --set-key-file <KEY> <FILE>      Set an internal key-value pair, with value being read from a file
      --split-channels[=<ENABLE>]      Split the grid such that each channel contains only a single PID combination [possible values: true, false]
      --upgrade[=<ENABLE>]             Convert the file format to the most recent version [possible values: true, false]
  -h, --help                           Print help
";

const CHANNEL_STR: &str = "c    entry        entry
-+------------+------------
0 1 × ( 2, -1) 1 × ( 4, -3)
1 1 × (21, -3) 1 × (21, -1)
2 1 × (22, -3) 1 × (22, -1)
3 1 × ( 2, 21) 1 × ( 4, 21)
4 1 × ( 2, 22) 1 × ( 4, 22)
";

const DEDUP_CHANNEL_DIFF_STR: &str = "b    x1               O(as^0 a^2)                       O(as^0 a^3)                       O(as^1 a^2)          
-+----+----+-----------+-----------+-------+-------------+-------------+-------+-----------+-----------+-------
0    2 2.25 6.5070305e2 6.5070305e2 0.000e0  -7.8692484e0  -7.8692484e0 0.000e0 1.1175729e2 1.1175729e2 0.000e0
1 2.25  2.5 5.9601236e2 5.9601236e2 0.000e0  -6.5623495e0  -6.5623495e0 0.000e0 1.0083341e2 1.0083341e2 0.000e0
2  2.5 2.75 5.1561247e2 5.1561247e2 0.000e0  -5.2348261e0  -5.2348261e0 0.000e0 8.9874343e1 8.9874343e1 0.000e0
3 2.75    3 4.1534629e2 4.1534629e2 0.000e0  -3.7590420e0  -3.7590420e0 0.000e0 7.3935106e1 7.3935106e1 0.000e0
4    3 3.25 3.0812719e2 3.0812719e2 0.000e0  -2.5871885e0  -2.5871885e0 0.000e0 5.6414554e1 5.6414554e1 0.000e0
5 3.25  3.5 2.0807482e2 2.0807482e2 0.000e0  -1.6762487e0  -1.6762487e0 0.000e0 3.9468336e1 3.9468336e1 0.000e0
6  3.5    4 9.6856769e1 9.6856769e1 0.000e0 -8.1027456e-1 -8.1027456e-1 0.000e0 1.9822014e1 1.9822014e1 0.000e0
7    4  4.5 2.2383492e1 2.2383492e1 0.000e0 -2.2022770e-1 -2.2022770e-1 0.000e0 5.3540011e0 5.3540011e0 0.000e0
";

const DEFAULT_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 7.5459110e2
1 2.25  2.5 6.9028342e2
2  2.5 2.75 6.0025198e2
3 2.75    3 4.8552235e2
4    3 3.25 3.6195456e2
5 3.25  3.5 2.4586691e2
6  3.5    4 1.1586851e2
7    4  4.5 2.7517266e1
";

const DELETE_BINS_02_57_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0 2.75    3 4.8552235e2
1    3 3.25 3.6195456e2
";

const DELETE_BINS_25_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 7.5459110e2
1 2.25  2.5 6.9028342e2
2  3.5    4 1.1586851e2
3    4  4.5 2.7517266e1
";

const DELETE_CHANNELS_STR: &str = "c    entry        entry
-+------------+------------
0 1 × ( 2, -1) 1 × ( 4, -3)
1 1 × (22, -3) 1 × (22, -1)
";

const KEY_VALUE_STR: &str = r"arxiv: 1505.07024
description: LHCb differential W-boson production cross section at 7 TeV
hepdata: 10.17182/hepdata.2114.v1/t4
initial_state_1: 2212
initial_state_2: 2212
key: value
lumi_id_types: pdg_mc_ids
mg5amc_repo: http://bazaar.launchpad.net/~maddevelopers/mg5amcnlo/3.1.2/
mg5amc_revno: 983
multiline: one
two
three
four
nnpdf_id: LHCBWZMU7TEV
pineappl_gitversion: v0.4.1-114-gdce19e0
results: ----------------------------------------------------------------------
   PineAPPL         MC        sigma      central         min      max 
                              1/100   sigma   1/1000   1/1000   1/1000
----------------------------------------------------------------------
 3.772955e+02  3.772821e+02   0.165   0.022   0.0357   0.0392   0.0313
 3.451417e+02  3.451342e+02   0.179   0.012   0.0217   0.0251   0.0172
 3.001260e+02  3.001231e+02   0.029   0.033   0.0096   0.0104   0.0076
 2.427612e+02  2.427624e+02   0.024   0.021   0.0049   0.0046   0.0060
 1.809773e+02  1.809799e+02   0.023   0.062   0.0143   0.0134   0.0154
 1.229334e+02  1.229354e+02   0.028   0.056   0.0157   0.0120   0.0200
 1.158685e+02  1.158603e+02   0.029   0.245   0.0708   0.0859   0.0514
 2.751727e+01  2.749798e+01   0.074   0.944   0.7014   0.7554   0.6281

runcard_gitversion: 82de4ad
x1_label: etal
x1_label_tex: $\eta_{\bar{\ell}}$
x1_unit: 
y_label: dsig/detal
y_label_tex: $\frac{\mathrm{d}\sigma}{\mathrm{d}\eta_{\bar{\ell}}}$
y_unit: pb
";

const MERGE_BINS_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 7.5459110e2
1 2.25  2.5 6.9028342e2
2  2.5 2.75 6.0025198e2
3 2.75    3 4.8552235e2
4    3 3.25 3.6195456e2
5 3.25  3.5 2.4586691e2
6  3.5  4.5 7.1692887e1
";

const REMAP_STR: &str = "b etal  x2  x3  dsig/detal 
   []   []  []     [pb]    
-+--+--+-+-+-+-+-----------
0  0  1 0 2 1 2 3.7729555e1
1  0  1 0 2 2 3 3.4514171e1
2  0  1 0 2 3 4 3.0012599e1
3  0  1 0 2 4 5 2.4276118e1
4  0  1 2 4 1 2 1.8097728e1
5  1  2 0 2 8 9 1.2293345e1
6  1  2 2 4 3 4 1.1586851e1
7  1  2 2 4 4 5 2.7517266e0
";

const REMAP_NO_REMAPPER_STR: &str = "Error: grid does not have a remapper
";

const REWRITE_CHANNELS_CONVOLVE_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 7.5534392e2
1 2.25  2.5 6.9342538e2
2  2.5 2.75 6.0526279e2
3 2.75    3 4.9140786e2
4    3 3.25 3.6782869e2
5 3.25  3.5 2.5085041e2
6  3.5    4 1.1874486e2
7    4  4.5 2.8214633e1
";

const REWRITE_CHANNELS_STR: &str = "c              entry                            entry                       entry                       entry                        entry                  entry
-+--------------------------------+-------------------------------+-----------------------+--------------------------------+-----------------------+---------------------
0 0.0000128881 × ( 2, -5)          0.050940490000000005 × ( 2, -3) 0.9490461561 × ( 2, -1) 0.0017222500000000003 × ( 4, -5) 0.9473907556 × ( 4, -3) 0.05089536 × ( 4, -1)
1 0.0017351381000000003 × (-5, 21) 0.9983312456 × (-3, 21)         0.9999415161 × (-1, 21)                                                          
2 1 × (22, -3)                     1 × (22, -1)                                                                                                     
3 0.9999995342 × ( 2, 21)          1.0000083656 × ( 4, 21)                                                                                          
4 1 × ( 2, 22)                     1 × ( 4, 22)                                                                                                     
";

const SCALE_BY_BIN_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 7.5459110e2
1 2.25  2.5 1.3805668e3
2  2.5 2.75 1.8007559e3
3 2.75    3 1.9420894e3
4    3 3.25 1.8097728e3
5 3.25  3.5 1.4752015e3
6  3.5    4 8.1107956e2
7    4  4.5 2.2013813e2
";

const SCALE_BY_ORDER_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 4.3317419e2
1 2.25  2.5 3.9555841e2
2  2.5 2.75 3.4506316e2
3 2.75    3 2.7972873e2
4    3 3.25 2.0918456e2
5 3.25  3.5 1.4266762e2
6  3.5    4 6.7845261e1
7    4  4.5 1.6435633e1
";

const SPLIT_CHANNELS_STR: &str = "c    entry
-+------------
0 1 × ( 2, -1)
1 1 × ( 4, -3)
2 1 × (21, -3)
3 1 × (21, -1)
4 1 × (22, -3)
5 1 × (22, -1)
6 1 × ( 2, 21)
7 1 × ( 4, 21)
8 1 × ( 2, 22)
9 1 × ( 4, 22)
";

const MULTIPLE_ARGUMENTS_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2  2.5 7.5454524e2
1  2.5    3 5.6456123e2
2    3 3.25 3.7400170e2
3 3.25  3.5 2.5300890e2
4  3.5    4 1.1909464e2
5    4  4.5 2.9004607e1
";

const ROTATE_PID_BASIS_NO_DIFF_STR: &str = "b    x1               O(as^0 a^2)                       O(as^0 a^3)                       O(as^1 a^2)          
-+----+----+-----------+-----------+-------+-------------+-------------+-------+-----------+-----------+-------
0    2 2.25 6.5070305e2 6.5070305e2 0.000e0  -7.8692484e0  -7.8692484e0 0.000e0 1.1175729e2 1.1175729e2 0.000e0
1 2.25  2.5 5.9601236e2 5.9601236e2 0.000e0  -6.5623495e0  -6.5623495e0 0.000e0 1.0083341e2 1.0083341e2 0.000e0
2  2.5 2.75 5.1561247e2 5.1561247e2 0.000e0  -5.2348261e0  -5.2348261e0 0.000e0 8.9874343e1 8.9874343e1 0.000e0
3 2.75    3 4.1534629e2 4.1534629e2 0.000e0  -3.7590420e0  -3.7590420e0 0.000e0 7.3935106e1 7.3935106e1 0.000e0
4    3 3.25 3.0812719e2 3.0812719e2 0.000e0  -2.5871885e0  -2.5871885e0 0.000e0 5.6414554e1 5.6414554e1 0.000e0
5 3.25  3.5 2.0807482e2 2.0807482e2 0.000e0  -1.6762487e0  -1.6762487e0 0.000e0 3.9468336e1 3.9468336e1 0.000e0
6  3.5    4 9.6856769e1 9.6856769e1 0.000e0 -8.1027456e-1 -8.1027456e-1 0.000e0 1.9822014e1 1.9822014e1 0.000e0
7    4  4.5 2.2383492e1 2.2383492e1 0.000e0 -2.2022770e-1 -2.2022770e-1 0.000e0 5.3540011e0 5.3540011e0 0.000e0
";

const ROTATE_PID_BASIS_DIFF_STR: &str = "b    x1                O(as^0 a^2)                          O(as^0 a^3)                          O(as^1 a^2)            
-+----+----+-----------+-----------+----------+-------------+-------------+----------+-----------+-----------+----------
0    2 2.25 6.5070305e2 6.5070305e2 -2.220e-16  -7.8692484e0  -7.8692484e0 -4.441e-16 1.1175729e2 1.1175729e2 -1.221e-15
1 2.25  2.5 5.9601236e2 5.9601236e2 -7.772e-16  -6.5623495e0  -6.5623495e0 -2.220e-16 1.0083341e2 1.0083341e2 -5.551e-16
2  2.5 2.75 5.1561247e2 5.1561247e2 -8.882e-16  -5.2348261e0  -5.2348261e0 -6.661e-16 8.9874343e1 8.9874343e1 -1.221e-15
3 2.75    3 4.1534629e2 4.1534629e2 -4.441e-16  -3.7590420e0  -3.7590420e0 -5.551e-16 7.3935106e1 7.3935106e1 -1.554e-15
4    3 3.25 3.0812719e2 3.0812719e2 -3.331e-16  -2.5871885e0  -2.5871885e0 -5.551e-16 5.6414554e1 5.6414554e1 -2.220e-16
5 3.25  3.5 2.0807482e2 2.0807482e2 -6.661e-16  -1.6762487e0  -1.6762487e0 -1.110e-16 3.9468336e1 3.9468336e1 -3.331e-16
6  3.5    4 9.6856769e1 9.6856769e1 -3.331e-16 -8.1027456e-1 -8.1027456e-1 -1.110e-16 1.9822014e1 1.9822014e1 -1.110e-15
7    4  4.5 2.2383492e1 2.2383492e1 -4.441e-16 -2.2022770e-1 -2.2022770e-1 -5.551e-16 5.3540011e0 5.3540011e0 -3.331e-16
";

const ROTATE_PID_BASIS_READ_CHANNELS_STR: &str = " c                 entry
---+-----------------------------------
0   0.013888888888888888 × (100, 100)
1   -0.020833333333333332 × (100, 103)
2   -0.006944444444444444 × (100, 108)
3   0.006944444444444444 × (100, 115)
4   0.004166666666666667 × (100, 124)
5   0.0027777777777777775 × (100, 135)
6   -0.013888888888888888 × (100, 200)
7   0.020833333333333332 × (100, 203)
8   0.006944444444444444 × (100, 208)
9   -0.006944444444444444 × (100, 215)
10  -0.004166666666666667 × (100, 224)
11  -0.0027777777777777775 × (100, 235)
12  0.020833333333333332 × (103, 100)
13  -0.0625 × (103, 103)
14  0.020833333333333332 × (103, 108)
15  0.010416666666666666 × (103, 115)
16  0.00625 × (103, 124)
17  0.004166666666666667 × (103, 135)
18  -0.020833333333333332 × (103, 200)
19  0.0625 × (103, 203)
20  -0.020833333333333332 × (103, 208)
21  -0.010416666666666666 × (103, 215)
22  -0.00625 × (103, 224)
23  -0.004166666666666667 × (103, 235)
24  0.006944444444444444 × (108, 100)
25  -0.020833333333333332 × (108, 103)
26  0.006944444444444444 × (108, 108)
27  0.003472222222222222 × (108, 115)
28  0.0020833333333333333 × (108, 124)
29  0.0013888888888888887 × (108, 135)
30  -0.006944444444444444 × (108, 200)
31  0.020833333333333332 × (108, 203)
32  -0.006944444444444444 × (108, 208)
33  -0.003472222222222222 × (108, 215)
34  -0.0020833333333333333 × (108, 224)
35  -0.0013888888888888887 × (108, 235)
36  -0.006944444444444444 × (115, 100)
37  -0.010416666666666666 × (115, 103)
38  0.024305555555555552 × (115, 108)
39  -0.003472222222222222 × (115, 115)
40  -0.0020833333333333337 × (115, 124)
41  -0.001388888888888889 × (115, 135)
42  0.006944444444444444 × (115, 200)
43  0.010416666666666666 × (115, 203)
44  -0.024305555555555552 × (115, 208)
45  0.003472222222222222 × (115, 215)
46  0.0020833333333333337 × (115, 224)
47  0.001388888888888889 × (115, 235)
48  -0.00625 × (124, 103)
49  -0.0020833333333333333 × (124, 108)
50  0.0020833333333333333 × (124, 115)
51  0.0012500000000000002 × (124, 124)
52  0.0008333333333333334 × (124, 135)
53  -0.004166666666666667 × (124, 200)
54  0.00625 × (124, 203)
55  0.0020833333333333333 × (124, 208)
56  -0.0020833333333333333 × (124, 215)
57  -0.0012500000000000002 × (124, 224)
58  -0.0008333333333333334 × (124, 235)
59  -0.004166666666666667 × (135, 103)
60  -0.0013888888888888887 × (135, 108)
61  0.0013888888888888887 × (135, 115)
62  0.0005555555555555556 × (135, 135)
63  -0.0027777777777777775 × (135, 200)
64  0.004166666666666667 × (135, 203)
65  0.0013888888888888887 × (135, 208)
66  -0.0013888888888888887 × (135, 215)
67  -0.0008333333333333334 × (135, 224)
68  -0.0005555555555555556 × (135, 235)
69  0.013888888888888888 × (200, 100)
70  0.004166666666666667 × (200, 124)
71  0.0027777777777777775 × (200, 135)
72  -0.013888888888888888 × (200, 200)
73  0.020833333333333332 × (200, 203)
74  0.006944444444444444 × (200, 208)
75  -0.006944444444444444 × (200, 215)
76  -0.004166666666666667 × (200, 224)
77  -0.0027777777777777775 × (200, 235)
78  -0.0625 × (203, 103)
79  -0.020833333333333332 × (203, 200)
80  0.0625 × (203, 203)
81  -0.020833333333333332 × (203, 208)
82  -0.010416666666666666 × (203, 215)
83  -0.00625 × (203, 224)
84  -0.004166666666666667 × (203, 235)
85  0.006944444444444444 × (208, 108)
86  0.003472222222222222 × (208, 115)
87  -0.006944444444444444 × (208, 200)
88  0.020833333333333332 × (208, 203)
89  -0.006944444444444444 × (208, 208)
90  -0.003472222222222222 × (208, 215)
91  -0.0020833333333333333 × (208, 224)
92  -0.0013888888888888887 × (208, 235)
93  0.024305555555555552 × (215, 108)
94  -0.003472222222222222 × (215, 115)
95  -0.0020833333333333337 × (215, 124)
96  -0.001388888888888889 × (215, 135)
97  0.006944444444444444 × (215, 200)
98  0.010416666666666666 × (215, 203)
99  -0.024305555555555552 × (215, 208)
100 0.003472222222222222 × (215, 215)
101 0.0020833333333333337 × (215, 224)
102 0.001388888888888889 × (215, 235)
103 0.004166666666666667 × (224, 100)
104 0.0020833333333333333 × (224, 115)
105 0.0012500000000000002 × (224, 124)
106 0.0008333333333333334 × (224, 135)
107 0.00625 × (224, 203)
108 0.0020833333333333333 × (224, 208)
109 -0.0020833333333333333 × (224, 215)
110 -0.0012500000000000002 × (224, 224)
111 -0.0008333333333333334 × (224, 235)
112 0.0027777777777777775 × (235, 100)
113 0.0013888888888888887 × (235, 115)
114 0.0008333333333333334 × (235, 124)
115 0.0005555555555555556 × (235, 135)
116 0.004166666666666667 × (235, 203)
117 0.0013888888888888887 × (235, 208)
118 -0.0013888888888888887 × (235, 215)
119 -0.0005555555555555556 × (235, 235)
120 0.16666666666666666 × (21, 100)
121 -0.25 × (21, 103)
122 -0.08333333333333333 × (21, 108)
123 0.08333333333333333 × (21, 115)
124 0.05 × (21, 124)
125 0.03333333333333333 × (21, 135)
126 -0.16666666666666666 × (21, 200)
127 0.25 × (21, 203)
128 0.08333333333333333 × (21, 208)
129 -0.08333333333333333 × (21, 215)
130 -0.05 × (21, 224)
131 -0.03333333333333333 × (21, 235)
132 0.16666666666666666 × (22, 100)
133 -0.25 × (22, 103)
134 -0.08333333333333333 × (22, 108)
135 0.08333333333333333 × (22, 115)
136 0.05 × (22, 124)
137 0.03333333333333333 × (22, 135)
138 -0.16666666666666666 × (22, 200)
139 0.25 × (22, 203)
140 0.08333333333333333 × (22, 208)
141 -0.08333333333333333 × (22, 215)
142 -0.05 × (22, 224)
143 -0.03333333333333333 × (22, 235)
144 0.25 × (103, 21)
145 0.08333333333333333 × (108, 21)
146 -0.08333333333333334 × (115, 21)
147 0.16666666666666666 × (200, 21)
148 -0.08333333333333334 × (215, 21)
149 0.05 × (224, 21)
150 0.03333333333333333 × (235, 21)
151 0.25 × (103, 22)
152 0.08333333333333333 × (108, 22)
153 -0.08333333333333334 × (115, 22)
154 0.16666666666666666 × (200, 22)
155 -0.08333333333333334 × (215, 22)
156 0.05 × (224, 22)
157 0.03333333333333333 × (235, 22)
";

const REWRITE_ORDER_CONVOLVE_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 1.8216658e2
1 2.25  2.5 1.6597039e2
2  2.5 2.75 1.4666687e2
3 2.75    3 1.2014156e2
4    3 3.25 9.0894574e1
5 3.25  3.5 6.2823156e1
6  3.5    4 3.0663454e1
7    4  4.5 7.8264717e0
";

const REWRITE_ORDER_READ_STR: &str = "o      order
-+----------------
0 O(as^1 a^1)
1 O(as^1 a^2)
2 O(as^1 a^2 lf^1)
3 O(a^3)
4 O(a^3 lf^1)
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["write", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn cc1() {
    let output = NamedTempFile::new("cc1.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--cc1",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DEFAULT_STR);
}

#[test]
fn cc2() {
    let output = NamedTempFile::new("cc2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--cc2",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DEFAULT_STR);
}

#[test]
fn delete_bins_02_57() {
    let output = NamedTempFile::new("deleted.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--delete-bins=0-2,5-7",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DELETE_BINS_02_57_STR);
}

#[test]
fn delete_bins_25() {
    let output = NamedTempFile::new("deleted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--delete-bins=2-5",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DELETE_BINS_25_STR);
}

#[test]
fn delete_channels() {
    let output = NamedTempFile::new("deleted3.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--delete-channels=1,3-4",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["read", "--channels", output.path().to_str().unwrap()])
        .assert()
        .success()
        .stdout(DELETE_CHANNELS_STR);
}

#[test]
fn key_value() {
    let output = NamedTempFile::new("set.pineappl.lz4").unwrap();
    let file = NamedTempFile::new("file").unwrap();

    file.write_str("one\ntwo\nthree\nfour").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--delete-key=runcard",
            "--set-key-value",
            "key",
            "value",
            "--set-key-file",
            "multiline",
            file.path().to_str().unwrap(),
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["read", "--show", output.path().to_str().unwrap()])
        .assert()
        .success()
        .stdout(KEY_VALUE_STR);
}

#[test]
fn merge_bins() {
    let output = NamedTempFile::new("bins.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--merge-bins=6-7",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(MERGE_BINS_STR);
}

#[test]
fn optimize() {
    // use `.pineappl` extension without `.lz4` to test `Grid::write` without compresssion
    let output = NamedTempFile::new("optimized.pineappl").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--optimize",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");
}

#[test]
fn remap() {
    let output = NamedTempFile::new("remapped.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--remap=0,1,2;0,2,4;1,2,3,4,5|:3|5:1,2,3,4,5,8,9|2:2",
            "--remap-norm-ignore=1",
            "--remap-norm=5",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(REMAP_STR);
}

#[test]
fn remap_norm_no_remapper() {
    let output = NamedTempFile::new("remapped.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--remap-norm=1",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(REMAP_NO_REMAPPER_STR);
}

#[test]
fn remap_norm_ignore_no_remapper() {
    let output = NamedTempFile::new("remapped.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--remap-norm-ignore=0",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(REMAP_NO_REMAPPER_STR);
}

#[test]
fn scale_by_bin() {
    let output = NamedTempFile::new("scale_by_bin.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--scale-by-bin=1,2,3,4,5,6,7,8",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(SCALE_BY_BIN_STR);
}

#[test]
fn scale_by_order() {
    let output = NamedTempFile::new("merged.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--scale-by-order=2,1,0.5,0.5",
            "--scale=0.5",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(SCALE_BY_ORDER_STR);
}

#[test]
fn split_channels() {
    let output = NamedTempFile::new("split-channels.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--split-channels",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DEFAULT_STR);

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["read", "--channels", output.path().to_str().unwrap()])
        .assert()
        .success()
        .stdout(SPLIT_CHANNELS_STR);
}

#[test]
fn dedup_channels() {
    let output = NamedTempFile::new("dedup-channels.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--split-channels",
            "--dedup-channels",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "diff",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DEDUP_CHANNEL_DIFF_STR);

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["read", "--channels", output.path().to_str().unwrap()])
        .assert()
        .success()
        .stdout(CHANNEL_STR);
}

#[test]
fn upgrade() {
    let output = NamedTempFile::new("upgraded.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--upgrade",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");
}

#[test]
fn multiple_arguments() {
    let output = NamedTempFile::new("multiple.merge.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "--merge-bins=0-1,2-3",
            "--scale=2",
            "--merge-bins=0-0",
            "--scale=0.5",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            output.path().to_str().unwrap(),
            "NNPDF40_nnlo_as_01180",
        ])
        .assert()
        .success()
        .stdout(MULTIPLE_ARGUMENTS_STR);
}

#[test]
fn rewrite_channels() {
    let output = NamedTempFile::new("ckm_channels.pineappl.lz4").unwrap();

    // 0 1 × ( 2, -1) 1 × ( 4, -3)
    // 1 1 × (21, -3) 1 × (21, -1)
    // 2 1 × (22, -3) 1 × (22, -1)
    // 3 1 × ( 2, 21) 1 × ( 4, 21)
    // 4 1 × ( 2, 22) 1 × ( 4, 22)

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--rewrite-channel", "0", "0.9490461561 * ( 2, -1) + 0.050940490000000005 * (2, -3) + 0.0000128881 * (2, -5) + 0.05089536 * (4, -1) + 0.9473907556 * (4, -3) + 0.0017222500000000003 * (4, -5)",
            "--rewrite-channel", "1", "0.9999415161 * (-1, 21) + 0.9983312456 * (-3, 21) + 0.0017351381000000003 * (-5, 21)",
            "--rewrite-channel", "3", "0.9999995342 * ( 2, 21) + 1.0000083656 * ( 4, 21)",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["read", "--channels", output.path().to_str().unwrap()])
        .assert()
        .success()
        .stdout(REWRITE_CHANNELS_STR);

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(REWRITE_CHANNELS_CONVOLVE_STR);
}

#[test]
fn rotate_pid_basis() {
    let pdg_to_pdg = NamedTempFile::new("pdg-to-pdg.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--rotate-pid-basis=PDG",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            pdg_to_pdg.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "diff",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            pdg_to_pdg.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(ROTATE_PID_BASIS_NO_DIFF_STR);

    let pdg_to_evol = NamedTempFile::new("pdg-to-evol.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--rotate-pid-basis=EVOL",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            pdg_to_evol.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "diff",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            pdg_to_evol.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
            "--ignore-channels",
        ])
        .assert()
        .success()
        .stdout(ROTATE_PID_BASIS_DIFF_STR);

    let evol_to_evol = NamedTempFile::new("evol-to-evol.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--rotate-pid-basis=EVOL",
            pdg_to_evol.path().to_str().unwrap(),
            evol_to_evol.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "diff",
            pdg_to_evol.path().to_str().unwrap(),
            evol_to_evol.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(ROTATE_PID_BASIS_NO_DIFF_STR);

    let evol_to_pdg = NamedTempFile::new("evol-to-pdg.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--rotate-pid-basis=PDG",
            // fix factors that are almost '1' to exact '1's
            "--rewrite-channel",
            "0",
            "1 * ( 2, -1) + 1 * ( 4, -3)",
            pdg_to_evol.path().to_str().unwrap(),
            evol_to_pdg.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "diff",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            evol_to_pdg.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(ROTATE_PID_BASIS_NO_DIFF_STR);

    let evol_to_evol_optimize = NamedTempFile::new("evol-to-evol-optimize.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--rotate-pid-basis=EVOL",
            // use the old name instead of `--split-channels` to test the alias
            "--split-lumi",
            "--optimize",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            evol_to_evol_optimize.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "read",
            // use the old name instead of `--channels` to test the alias
            "--channels",
            evol_to_evol_optimize.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout(ROTATE_PID_BASIS_READ_CHANNELS_STR);
}

#[test]
fn rewrite_order() {
    let output = NamedTempFile::new("rewrite-order.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--rewrite-order",
            "0",
            "as1a1",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["read", "--orders", output.path().to_str().unwrap()])
        .assert()
        .success()
        .stdout(REWRITE_ORDER_READ_STR);

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(REWRITE_ORDER_CONVOLVE_STR);
}
