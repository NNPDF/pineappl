# Tutorial for PineAPPL's CLI

Welcome to PineAPPL's CLI tutorial! Here we explain the basics of PineAPPL's
command-line interface (CLI): that's the program `pineappl` that you can you use
inside your shell. This will also introduce and explain the terminology needed
to understand the C, Python and Rust API.

First, let's explain what interpolation grids are and what PineAPPL is. If you
already know that, feel free to skip these sections!

### What are interpolation grids?

Interpolation grids store theoretical predictions that are not convoluted with
PDFs yet. After they've been generated once, the grids can be convoluted with
arbitrary PDFs in a fraction of a second. This is very advantageous for the
following applications:

- *PDF set dependence* of a prediction. Are the predictions with CT18, MSHT20
  and NNPDF3.1 compatible? If they are not, which part(s) of their PDF sets is
  responsible for this?
- *PDF fits* require the evaluation and fitting of theory predictions to
  experimental data with changing PDFs, which are impossible without
  PDF-independent theory predictions.

### What is PineAPPL?

PineAPPL is one interpolation library, and provides interfaces to numerical
calcuations to generate interpolation grids. It is not the only library that
does that, but in fact there at least two others: [APPLgrid] and [fastNLO]. A
distinguishing feature of PineAPPL is

- support for arbitrary fixed-order calculations, in particular EW and mixed
  QCD-EW corrections. Furthermore it aims to
- provide *excellent* tooling, which is easy to use and provides maximum
  insight into theory predictions.

## Installation

To take the tutorial, you'll need PineAPPL's CLI; follow the installation
section in the [README](../README.md). Next, we'll need a fresh directory that
we can work in, for instance,

    cd $(mktemp -d)

will create a temporary directory. Finally, we'll need a grid,

    wget 'https://github.com/N3PDF/pineappl/raw/master/pineappl_cli/data/LHCB_WP_7TEV.pineappl.lz4'

which we'll use together with the CLI.

## Performing convolutions

Now that we have a grid, let's perform a convolution with a PDF set:

    pineappl convolute LHCB_WP_7TEV.pineappl.lz4 CT18NNLO

We've chosen the default CT18 PDF set to use in this tutorial, because it's the
easiest to type. If you get this error,

    error: Invalid value for '<PDFSETS>...': The PDF set `CT18NNLO` was not found

install the PDF set in LHAPDF, or use a different PDF set—the numbers won't
matter for the sake of the tutorial. If the command was successful, you should
see the following output:

    LHAPDF 6.4.0 loading /home/cschwan/prefix/share/LHAPDF/CT18NNLO/CT18NNLO_0000.dat
    CT18NNLO PDF set, member #0, version 1; LHAPDF ID = 14000
    bin   etal    disg/detal  scale uncertainty
    ---+----+----+-----------+--------+--------
      0    2 2.25 3.8562616e2   -3.69%    2.65%
      1 2.25  2.5 3.5553189e2   -3.70%    2.73%
      2  2.5 2.75 3.0928086e2   -3.68%    2.79%
      3 2.75    3 2.4991572e2   -3.66%    2.83%
      4    3 3.25 1.8610527e2   -3.62%    2.85%
      5 3.25  3.5 1.2614289e2   -3.57%    2.87%
      6  3.5    4 5.9206239e1   -3.47%    2.84%
      7    4  4.5 1.3881523e1   -3.29%    2.75%
    Thanks for using LHAPDF 6.4.0. Please make sure to cite the paper:
      Eur.Phys.J. C75 (2015) 3, 132  (http://arxiv.org/abs/1412.7420)

On your computer the output will be slightly different depending on your LHAPDF
installation. If you don't want to see LHAPDF messages, add the option
`--silence-lhapdf` after `pineappl` and before `convolute`:

    pineappl --silence-lhapdf convolute LHCB_WP_7TEV.pineappl.lz4 CT18NNLO

Let's have a closer look at what the output shows us:

- the `bin` column shows that there are 8 bins, labeled 0 to 7,
- the next two columns `etal` shows the observable used to define the bins and
  its left and right bin limits,
- `disg/detal` has a spelling mistake unfortunately and should read
  `dsig/detal`. It shows the differential cross section (`dsig`) for the bin
  limits given in the previous two columns. Finally,
- the last two columns show the scale uncertainty, which typically is
  asymmetric.

## Getting help

One of the most difficult aspects of learning a new program is remembering
*how* to achieve certain tasks and *what* to type. Fortunately, this is easy
with `pineappl`. If you don't want to remember a bunch of commands just
remember that you can type

    pineappl

which is a shortcut for

    pineappl --help

This will give you a list of program options and *subcommands*, under which all
operations in `pineappl` are grouped, and a corresponding description. You'll
be familiar with the concept of subcommand if you're using `git`: `add`,
`commit` and `push` are well-known subcommands of it.

To get more help on a specific subcommand, for instance `convolute`, which
we've used already, run

    pineappl convolute --help

Depending on the version of PineAPPL this will show an output similar to the
following:

    pineappl-convolute
    Convolutes a PineAPPL grid with a PDF set

    USAGE:
        pineappl convolute [OPTIONS] <INPUT> <PDFSETS>...

    ARGS:
        <INPUT>         Path of the input grid
        <PDFSETS>...    LHAPDF id(s) or name of the PDF set(s)

    OPTIONS:
        -a, --absolute              Show absolute numbers of the scale variation
        -b, --bins <BINS>...        Selects a subset of bins
        -h, --help                  Print help information
        -i, --integrated            Show integrated numbers (without bin widths) instead of differential
                                    ones
        -o, --orders <ORDERS>...    Select orders manually
        -s, --scales <SCALES>       Set the number of scale variations [default: 7] [possible values: 1,
                                    3, 7, 9]

If you read the help message carefully, you'll for instance notice that the
scale uncertainty shown previously is a 7-point variation, because the default
value of `--scales` is `7`.

## What does this grid contain?

If you're experienced, you've already inferred from the file name of the grid
and the observable name what the convoluted numbers will most likely show.
However, how do we make certain? Specifically, we'd like to know the answers to
the following questions:

1) for which process is the prediction for? Which observable is shown?
2) Is there a corresponding measurement? Where's the paper for that
   measurement? Where's can we get the measurements corresponding to the
   predictions?
3) Which program/generator was used to produce this grid? What version? Which
   runcards were used, which parameter values?

These questions do not concern the main data, the theory predictions, but
instead the *metadata*, that's the data about the data. The subcommand that
we'll need is `info`. The easiest way to answer all the questions is to print
all metadata using

    pineappl info --show LHCB_WP_7TEV.pineappl.lz4

This will print out very long list of key–value pairs, of which the following
are important:

- `description` gives a short description of the process.
- The value of the key `arxiv` gives us the corresponding paper, in this case
  for the experimental measurement.
- The measured data itself we can get from the location specified by `hepdata`,
  which is a digital object identifier (DOI). Go to <https://www.doi.org/> and
  enter the string there. This will get you to the right table on
  <https://www.hepdata.net> for the corresponding observable.
- the presence of `mg5amc_repo` and `mg5amc_revno` signal that
  [Madgraph5_aMC@NLO] was used to calculate the predictions, and `runcard` are
  the runcards used to run the calculation.

[APPLgrid]: https://applgrid.hepforge.org/
[fastNLO]: https://fastnlo.hepforge.org/
[Madgraph5_aMC@NLO]: https://launchpad.net/mg5amcnlo
