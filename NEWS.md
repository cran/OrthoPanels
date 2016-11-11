# OrthoPanels 1.1.0

## New features

OrthoPanels can now handle cases that enter the panel late
(refreshment), not just drop outs. Late entrants can either have NAs
for time entries prior to their joining, or not be present in the data
at early times (possible only when using the formula interface).

Add a dataset of the dynamics of labour demand in the United Kingdom,
based on Arellano and Bond (1991)

## Bug fixes and minor improvements

- missing data only in the first wave of `X` doesn't affect the result.

- handle case and time variables that aren't `1:N` and `1:T`
