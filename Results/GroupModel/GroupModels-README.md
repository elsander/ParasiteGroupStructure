This folder contains results used in the manuscript, and were produced
by the code in `../../Code/GroupModelAlgorithm/`.

## Folder Names

The folder names are structured with the following information,
separated by dashes:
- food web name
- flag for concomitant predation
- flag for degree correction

## File Names

File names in each folder are structure with the following
information, separated by dashes:
- file name of food web from `../../Data/`
- "G-" followed by the maximum number of groups
- "DC-" followed by a flag for degree correction
- "alpha-" followed by the alpha parameter used for the gamma prior
(NA for the uncorrected model)
- "beta-" followed by the beta parameter used for the gamma prior
(NA for the uncorrected model)
- "Marginal" followed by the marginal likelihood of the grouping (here, the
  dash results from the fact that the marginal likelihood is negative)

## File Structure

Each file contains a row of *N* space-separated integers, where *N* is
the number of species in the food web. Each integer represents group
identity for a species, ordered as they are ordered in the
corresponding food web matrix.
