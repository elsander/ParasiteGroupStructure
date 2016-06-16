# Calculating the imbalance of ecological categories across a given partitioning
Author: Matt Michalska-Smith
Requirements: julia 0.4.5, Iterators package

The file `Imbalance.jl` contains functions for the analysis, while the file
`RunImbalance.jl` calls these functions, checks for errors, and writes the
results to file.

In `Imbalance.jl`:

| Function | Description | Inputs | Outputs |
|----------|----------------------|------------------------|--------|
| `Imbalance` | Calculates the imbalance of a given ecological category distribution | a presence vector array* | the (single) imbalance value |
| `Pvalue` | Calculates the probability of obtaining a given ecological category distribution | a presence vector array | the (single) p-value |
| `GetPotentialPresences` | Prunes all possible ecological category distributions across a given set of groups for feasibility | vector of group sizes and the (integer) number of species in the focal ecological category | a presence vector array |
| `GetCompatiblePresences` | Prunes all possible combinations of ecological category distributions for feasibility | vector of group sizes and a presence vector array | a presence vector array |
| `GetAllImbalances` | Constructs a table indicating the pvalue of each possible imbalance value | a presence vector array | a table of imbalances and p-values |
| `Emppvalue` | Calculates the pvalue of the empirically observed distribution (sums p-values for greater imbalances) | an empirical distribution and a table of imbalances and p-values | a single p-value |
| `SampleDists` | Samples several ecological category distributions for a given partitioning and calculates their imbalances to estimate the empirical pvalue | a focal ecological category, the empirical imbalance value, a partitioning vector, a vector of ecological categories, the number of samples to run, and the name of the file the results should be written to | None |

To run the analysis for a given ecological categorization and partitioning, run

    julia RunImbalance.jl [DATAFILE] [numsamps] [FOCAL_CATEGORY] [OUTDIR]

from the command line. Where `DATAFILE` is a file with two columns (separated by
whitespace): the first being the ecological categorization and the second the
partitioning; `NUMSAMPS` is the (positive integer) number of samples to run or
"-1" to indicate analytic calculation of the p-value; `FOCAL_CATEGORY` is
either one of the values in the first column of `DATAFILE` or "all" (indicating)
that an overall imbalance value should be calculated, rather than a pairwise one;
and `OUTDIR` is the directory in which to put all output files.
