This folder contains results used in the manuscript, and were produced
by the code in `../../Code/Imbalance/`.

## File Names

The file name is structured with the following information, separated
by dashes:
- food web name
- flag for concomitant predation
- flag for degree correction
- number of groups
- ecological category

Each file corresponds to a grouping in `../GroupModel/`.

## File Structure

Each file has four columns:
- ecological category
- imbalance value
- estimated p-value
- number of samples

Files have 100 lines each, showing the estimated p-value as the number
of samples increases.
