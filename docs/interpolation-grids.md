# Interpolation grids and PineAPPL

This document explains what interpolation grids are and what PineAPPL is.

## What are interpolation grids?

Interpolation grids store theoretical predictions that are not convoluted with
PDFs yet. After the predictions have been generated once, the grids can be
convoluted with arbitrary PDFs in a fraction of a second. This is very
advantageous for the following applications:

- *PDF set dependence* of a prediction. Are the predictions with CT18, MSHT20
  and NNPDF3.1 compatible? If they are not, which parts of their PDF sets are
  responsible for this?
- *PDF fits* require the evaluation and fit of theory predictions to
  experimental data with changing PDFs, which are impossible without
  PDF-independent theory predictions.

## What is PineAPPL?

PineAPPL is an interpolation library, and provides interfaces to numerical
calculations to generate interpolation grids. It is not the only library that
does that, as there at least two others: [APPLgrid] and [fastNLO].
Distinguishing features of PineAPPL are

- support for arbitrary fixed-order calculations, in particular EW and mixed
  QCD–EW corrections. Furthermore it aims to
- provide *excellent* tooling, which is easy to use and provides maximum
  insight into theory predictions.

[APPLgrid]: https://applgrid.hepforge.org/
[fastNLO]: https://fastnlo.hepforge.org/
