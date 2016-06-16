#ParasiteGroupStructure: code and data to support and replicate Sander
# *et al.* (*in prep.*)

This repository contains the code, data, and results needed to
replicate the findings in Sander *et al.* (*in prep.*).

The `Data/` folder contains data from Dunne *et al.* 2013 (PLoS
Biology), formatted for use in Sander *et al.*. The `Results/` folder
contains the results used in analyses (groupings, imbalance values,
and imbalance p-values). The `Code/` contains programs and scripts
used in our analyses: the group model search algorithm, imbalance and
p-value calculation, mutual information, and statistical
analyses. More detailed information about the files is given in README
files in repository subfolders. For a description of the methodology
and results, see Sander *et al.* (*in prep.*).

Note that throughout, flags are used to distinguish between
degree-corrected and uncorrected models. A 0 flag corresponds to the
uncorrected model, and a 1 flag corresponds to the degree-corrected
model. Similarly, there is a flag for whether or not concomitant
predation is included. "Par" is used for no concomitant predation, and
"ParCon" for concomitant predation.
