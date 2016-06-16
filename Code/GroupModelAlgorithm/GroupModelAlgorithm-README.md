Run multiple-chain Markov Chain Monte Carlo to partition nodes
of a network into a set number of groups, with the option to
correct for degree

Author: Elizabeth Sander (esander(at)uchicago(dot)edu)
Version: 1.0.0

Requirements: gcc compiler, GSL and GSL CBLAS libraries for C

Inputs:
- number of species
- file path to an adjacency matrix
- random (integer) seed
- number of steps for each chain to take
- number of chains
- maximum number of groups
- 0/1 flag for degree correction
- alpha parameter for gamma prior for the rate parameter omega
(degree corrected model only)
- beta parameter for gamma prior for the rate parameter omega
(degree corrected model only)

Example function call from terminal, if NOT correcting for degree:
./FindGroups N ./path/to/file seed steps nchains maxgroups degreeflag
Example function call from terminal, if correcting for degree:
./FindGroups N ./path/to/file seed steps nchains maxgroups degreeflag alpha beta

A chain swap is attempted every 20 steps. This is hardcoded in the
global variables in SearchAlgs.c

Prints a vector of length N to file. This vector contains the best
partition found by the algorithm. Each integer i in the vector
represents the group identity of species i.
The output file will be named:
"<./path/to/file>-G-<maxgroups>-DC-<degreeflag>-alpha-<alpha>-beta-<beta>-Marginal<marginal loglikelihood>.txt"

Therefore, you should copy your adjacency matrix to the folder where
you want your output file.
Note that if not using the degree-corrected model, <alpha> and
<beta> in the output file name will be NA.
