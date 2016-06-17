Author: Liz Sander (esander at uchicago dot edu)
Requirements: R >= 3.3.0, dplyr, agricolae, stringr, vegan, tidyr

## Statistical Analysis of Results

The main results are in `StatisticalAnalysis.Rmd`. This file can be
compiled into markdown or html by
[knitting it](http://kbroman.org/knitr_knutshell/pages/Rmarkdown.html). Code
chunks can also be run separately in `R`.

`ImbalanceTableSampling.R` contains code to create a data frame of
imbalance values from the sampling-based imbalance results in a
folder.

`MutualInformation.R` contains code to calculate the mutual
information shared by two partitions.

`CalculateEvenness.R` contains code to calculate evennness for all
partitions in the `../../Results/GroupModel/` directory.

Note that the analysis scripts assume that only the best partition is
kept; that is, a folder will only contain one partition for each
number of groups.
