Metadata
========

Metadata is stored in PineAPPL grids in the form of key–value pairs, which can
be read indiviually with `pineappl read --get <KEY> <GRID>`, and `pineappl read
--show <GRID>` shows all key–value pairs.

Known keys
----------

- `arxiv`
- `description`
- `hepdata`
- `inspire`
- `results`
- `results_pdf`

CLI-recognized keys
-------------------

The following keys are used in the CLI when printing numbers resulting from
convolutions:

- `x1_label`: label of the first dimension for every bin. If the bins have more
  than one dimensions, keys with higher indices like `x2_label` and `x3_label`
  are used.
- `x1_unit`: the physical unit for the first dimension for every bin. If the
  bins have more than one dimension, keys with higher indices like `x2_unit`
  and `x3_unit` are used.
- `y_label`: label of the quantities stored in this grid.
- `y_unit`: physical unit of quantities stored in this grid.

For each missing key a default value will be used

Physical units
--------------

Recognized units are:

- `pb`, `fb`: picobarn, femtobarn
- `GeV`: gigaelectronvolt

`pineappl plot ...`-recognized keys
-----------------------------------

- `x1_label_tex`:
- `y_label_tex`
