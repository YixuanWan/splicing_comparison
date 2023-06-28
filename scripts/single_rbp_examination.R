library(dplyr)
library(tidyr)
library(as.data.table)
library(ggplot2)

rbp_deseq <- function(file){
  if("symbol" %in% colnames(file)){
    df = file |> 
      filter(symbol %in% RBP) |> 
      dplyr::select(log2fold_change, padj, deseq2_table_name, symbol)
  }
  else{
    df = file |> 
      filter(gene_name %in% RBP) |> 
      dplyr::select(log2fold_change, padj, deseq2_table_name, gene_name) 
  }
  return(df)
}

# examine KALRN PSI in in vitro experiments
splicingTable = read.csv("/Users/Ewann/splicing_comparison/data/majiq/splicing_full_delta_psi_tables.csv")
experiments = read.csv("/Users/Ewann/splicing_comparison/samplesheet/tdp_experiments_updated-tdp_experiments_updated.csv")

kalrn = c("chr3:124701255-124702038", "chr3:124701598-124702038", "chr3:124700033-124700975", "chr3:124700033-124701093")

splicingTable |> 
  filter(paste_into_igv_junction %in% kalrn) |> 
  select(gene_name,paste_into_igv_junction,junc_cat,baseline_PSI,contrast_PSI,comparison,strand, junc_cat) |>  
  mutate(is_cryptic = contrast_PSI > 0.1 & baseline_PSI < 0.05)  |> 
  filter(paste_into_igv_junction == "chr3:124701255-124702038") |> 
  left_join(experiments |> select(comparison, experiment)) |> 
  mutate(experiment = ifelse(is.na(experiment), comparison, experiment)) |> 
  mutate(experiment = forcats::fct_reorder(experiment, contrast_PSI)) |> 
  ggplot(aes(x = contrast_PSI, y = comparison, fill = is_cryptic)) +
  geom_col() +
  theme_minimal() + 
  labs(
    title = "KALRN_chr3:124701255-124702038",
    subtitle = "PSI in in vitro KD studies",
    x = "PSI",
    y = "Experiment",
    fill = "Cryptic"
  ) + 
  facet_wrap(~junc_cat)

# RBP cell-type differences check
filepath = list.files("/Users/Ewann/splicing_comparison/data/deseq2", full.names = TRUE)
suffix = ".csv"
filenames = list.files("/Users/Ewann/splicing_comparison/data/deseq2")
experiment_names = gsub(suffix,"",filenames)
estimate_files = purrr::map(filepath,my_clean_reader)
estimate_files = purrr::map2(estimate_files, experiment_names, ~cbind(.x, deseq2_table_name = .y))

RBP = c("SRSF3", "SNRNP70", "ABCF1", "KALRN")
deseq2 = purrr::map(estimate_files, rbp_deseq) |> 
  rbindlist(use.names = FALSE) 

fc_rbp = deseq2 |> 
  pivot_wider(names_from = "gene_name", values_from = log2fold_change, id_cols = deseq2_table_name, names_prefix = "log2fold_change_")

deseq2 = deseq2 |> 
  pivot_wider(names_from = "gene_name", values_from = padj, id_cols = deseq2_table_name, names_prefix = "padj_") |> 
  full_join(fc_rbp)

rbp_log2fc = experiments |> 
  left_join(deseq2) |> 
  select(comparison, cell.type, experiment, contains("log2fold_change")) |> 
  pivot_longer(cols = starts_with("log2"), names_to = "rbp", names_prefix = "log2fold_change_", values_to = "log2foldChange") 

rbp_deseq2 = experiments |> 
  left_join(deseq2) |> 
  select(comparison, experiment, cell.type, contains("padj")) |> 
  pivot_longer(cols = starts_with("padj"), names_to = "rbp", names_prefix = "padj_", values_to = "padj") |> 
  right_join(rbp_log2fc)

GENE = "tdp"
rbp_deseq2 |> 
  filter(rbp == GENE) |> 
  filter(!grepl("cycloheximide", comparison)) |> 
  filter(!grepl("upf1", comparison)) |> 
  mutate(comparison = forcats::fct_reorder(comparison, log2foldChange)) |> 
  ggplot(aes(x = comparison, y = log2foldChange, fill = cell.type)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title = GENE, 
    x = "Comparison",
    y = "log2FoldChange",
    fill = "Cell type"
  )
