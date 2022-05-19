# CLI reference

## `REMAPPING`: Remapping parameter specification

This section specifies the `REMAPPING` parameter of `pineappl remap`.

### Motivation

For performance/simplicity reasons most Monte Carlo generators (and also the
PineAPPL `Grid::fill` method) neither support

1) multi-dimensional distributions nor
2) distributions whose bin sizes are not equally sized *during generation*.

To work around this problem a grid with a one-dimensional distribution with
equally-sized bins can be generated instead, and afterwards the bins can be
'remapped' to an N-dimensional distribution using the limits specified with the
`REMAPPING` string.

### Reference

The remapping string uses the following special characters to achieve this
(note that some of these characters must be escaped in certain shells):

- `,`: The comma constructs 1-dimensional bin limits (1DBL). For example,
  the 1DBL `0,0.2,0.4,0.6,1` expects the grid to have 4 bins whose bin limits
  will be (0-0.2), (0.2-0.4), (0.4-0.6) and (0.6,1).
- `;`: If higher-dimensional bins are needed, the n-dimensional bin limits
  (NDBL) are constructed from a Cartesian product of 1DBL separated with a
  semicolon. For example, `0,0.5,1;0,1,2` expects the grid to have 4 bins,
  whose 2DBL will be are: (0-0.5;0-1), (0-0.5;1-2), (0.5-1;0-1) and
  (0.5-1;1-2).
- `|`: The previous operators are enough to construct NDBL with
  differently-sized bins, but they can not construct the following bin limits:
  (0-1;0-1), (0-1;1-2), (1-2;0-2), (1-2;2-4), (1-2;4-6); here the 1DBL for the
  second dimension depend on the first dimension and also have a different
  number of bins. For the first two bins the 1DBL is `0,1,2`, but for the last
  three bins the 1DBL are `0,2,4,6`. This can be achieved using the following
  remapping string: `0,1,2;0,1,2|0,2,4,6`. Note that there have to be two 1DBL
  separated by `|`, because the first dimension has two bins. If there are more
  dimensions and/or bins, the number of 1DBL separated by `|` must match this
  number accordingly. An example of this is the following remapping string:
  `0,1,2;-2,0,2;0,1,2|1,2,3|2,3,4|3,4,5|4,5,6|5,6,7`. Here the third dimension
  has 6 1DBL separated by `|` because the first dimension has 2 bins and the
  second dimension has 3 bins.

  If the 1DBL is an empty string, the previous 1DBL is repeated, for example
  `0,1,2;0,1,2;0,1,2||0,2,4` is shorthand for `0,1,2;0,1,2;0,1,2|0,1,2|0,2,4`.
- `:`: The last feature of `|` can combined with `:`, which is used to 'cut'
  out bins from the left and/or right. For example, the remapping string
  `0,1,2;0,1,2,3:2|:1||:1|2:` is a more succinct way of writing the following
  remapping string: `0,1,2;0,1|0,1,2|0,1,2,3|0,1,2|2,3`.

Finally note that the differential cross sections are calculated using the bin
sizes (the product of bin widths of each dimension) given by the remapping
string. The option `--ignore-obs-norm` can be used to remove certain dimensions
from the bin size determination, for example `'0,10,20;0,2,4' --ignore-obs-norm
1` will normalize the bins with a size of `2` because the first dimension (with
index `1` will be ignored)
