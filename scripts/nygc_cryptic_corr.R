library(arrow)
library(dplyr)
library(tidyr)
library(data.table)


# load in dataset
nygc = read_parquet("/Users/Ewann/splicing_comparison/data/selective_cryptic_psi_in_nygc.parquet")

psi_tdp = nygc |> 
  filter(tdp_path == "path") |> 
  select(psi, paste_into_igv_junction) |> 
  spread(key = paste_into_igv_junction, value = psi)
