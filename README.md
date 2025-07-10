# SNPhub

The goal of SNPhub is to …

## Installation

You can install the development version of SNPhub like so:

``` r
remotes::install_github("tbergmann92/SNPhub")
```

## Example
``` r
library(SNPhub)
# import data
mat <- read.csv("inst/extdata/RAW_Genotype_Data_Subset.csv", sep = ",", quote = "", row.names = 1)
mat <- as.matrix(mat)

# run snphub
snps <- snphub(mat)
```