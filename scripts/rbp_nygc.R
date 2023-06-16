library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

  

# import files
human_rbp = fread("/Users/Ewann/splicing_comparison/data/human_rbps.csv",header=TRUE)
tpm_nygc = fread("/Users/Ewann/splicing_comparison/data/nygc/rsem_tpm_nygc.csv", header = TRUE)

rbp_tpm = tpm_nygc |> 
  mutate(gene = gsub("\\..*", "", gene)) |> 
  right_join(human_rbp, by = c("gene" = "ensmbleID")) |> 
  select(-gene) |> 
  tibble::column_to_rownames(var = "gene_name") |> 
  t() |> as.data.table()



pdf("/Users/Ewann/splicing_comparison/results/nygc/rbp_tdp_corr.pdf")
for (r in colnames(rbp_tpm)){
  plot = rbp_tpm |> 
    ggplot(aes(x = TARDBP, y =rbp_tpm[[r]])) + 
    geom_point() +
    geom_smooth(method = "lm") + 
    ggpubr::stat_cor() +
    theme_minimal() +
    labs(
      title = "NYGC TPM",
      y = r)
  print(plot)
}
dev.off()

# single RBP expression examination
RBP = "SNRNP70"
ensembl_id = "ENSG00000104852"
single_rbp = tpm_nygc |> 
  mutate(gene = gsub("\\..*", "", gene)) |> 
  filter(gene == ensembl_id)|>
  mutate(gene = RBP) |> 
  tibble::column_to_rownames(var = "gene") |> 
  t() |> 
  as.data.table() |> 
  cbind(rbp_tpm) 

single_rbp |> 
  ggplot(aes(x = TARDBP, y = single_rbp[[RBP]])) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  ggpubr::stat_cor() +
  theme_minimal() +
  labs(
    title = "NYGC TPM",
    y = RBP)

