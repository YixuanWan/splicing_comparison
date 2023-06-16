library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

# read MAJIQ files
splicing_full = fread("/Users/Ewann/splicing_comparison/data/majiq/splicing_full_delta_psi_tables.csv")
# load in DESEQ2 files