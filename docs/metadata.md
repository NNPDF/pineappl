# Metadata

Metadata is stored in PineAPPL grids in the form of key--value pairs, which can
be read indiviually with `pineappl read --get <KEY> <GRID>`, and `pineappl read
--show <GRID>` shows all key--value pairs.

## Known keys

- `arxiv`
- `description`
- `hepdata`
- `inspire`
- `results`
- `results_pdf`

## CLI-recognized keys

The following keys are used in the CLI when printing numbers resulting from
convolutions:

- `x1_label`: label of the first dimension for every bin. If the bins have more
  than one dimension, keys with higher indices like `x2_label` and `x3_label`
  are used.
- `x1_unit`: the physical unit for the first dimension for every bin. If the
  bins have more than one dimension, keys with higher indices like `x2_unit`
  and `x3_unit` are used.
- `y_label`: label of the quantities stored in this grid.
- `y_unit`: physical unit(s) of the quantities stored in this grid. If
  differential cross sections are stored, also the unit of the dividing
  observable must be given, for instance a cross sections differential in an
  invariant mass could have the units `fb/GeV`.

In all values avoid long strings and try to avoid using spaces, for instance
use `d2sig/dy/dMll`

For each missing key a default value will be used.

### Recommended values



### Physical units

Recognized units are:

- `pb`, `fb`: picobarn, femtobarn
- `GeV`: gigaelectronvolt

Ratios are denoted using `/`, and powers are denoted using `^2`

## `pineappl plot ...`-recognized keys

These should contain the equivalent of `x1_label` with (La)TeX commands. If
symbols from math mode are used, they must be enclosed with `$`.

- `x1_label_tex`
- `y_label_tex`
