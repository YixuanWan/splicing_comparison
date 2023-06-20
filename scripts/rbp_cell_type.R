library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

# read in nygc data
cell.type <- fread("/Users/Ewann/splicing_comparison/data/nygc/All_1917_samples_Darmanis_dtangle_deconv.tsv")
tpm_nygc = fread("/Users/Ewann/splicing_comparison/data/nygc/rsem_tpm_nygc.csv", header = TRUE)

nygc_by_gene = tpm_nygc |> 
  mutate(gene = gsub("\\..*", "", gene)) |> 
  distinct(gene, .keep_all = TRUE) |> 
  tibble::column_to_rownames(var = "gene") |> 
  t() |> as.data.table(keep.rownames = TRUE) |> 
  rename(sample = rn)

metadata = cell.type |> 
  select(sample, cell, deconv) |> 
  unique() |> 
  pivot_wider(id_cols = "sample", names_from = "cell", values_from = "deconv") 
  

# add rbp and tdp columns
tmp = nygc_by_gene |> 
  select("ENSG00000104852", "ENSG00000112081", "ENSG00000204574", "ENSG00000160145", "ENSG00000120948", sample) 
colnames(tmp) = c("SNRNP70", "SRSF3", "ABCF1", "KALRN", "TARDBP", "sample")

metadata = metadata |> 
  left_join(tmp) |> 
  column_to_rownames(var = "sample")

# perform regression
plot(metadata)

multiple_regression <- lm(SRSF3 ~ TARDBP + astrocytes + microglia, data = metadata)
summary(multiple_regression)
