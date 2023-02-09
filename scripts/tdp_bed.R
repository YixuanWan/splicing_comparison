library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(data.table)
library(plyranges)

import_bed <- function(file){
  df = data.table::as.data.table(import.bed(file))
  return(df)
}



# import files
suffix = '.bed'
filepath = list.files("/Users/Ewann/splicing_comparison/data/tdp_bed",
                   pattern = suffix,
                   full.names = TRUE)

tdp_bed = purrr::map(filepath, import_bed)

# get sample id
id = purrr::map(filepath, basename)
id = gsub(suffix, "", id)

# fill <name>with sample id
tdp_bed = purrr::map2(tdp_bed, id, ~cbind(.x, name2 = .y))

# combine them into a long 
metatable = data.table::rbindlist(tdp_bed)

meta_bed = metatable |> 
  select(-name) |> 
  mutate(name = name2, .keep = "unused") 

# write out the file
rtracklayer::export(meta_bed, "/Users/Ewann/splicing_comparison/data/tdp_bed/TARDBP_all.bed")


# to merge overlapping ranges
all_bed = import.bed("/Users/Ewann/splicing_comparison/data/tdp_bed/TARDBP_all.bed")
reduced = reduce(all_bed)
rtracklayer::export(reduced, "/Users/Ewann/splicing_comparison/data/tdp_bed/TARDBP_reduced.bed")

