library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

# load in files
tpm_nygc = fread("/Users/Ewann/splicing_comparison/data/nygc/rsem_tpm_nygc.csv", header = TRUE)
expression_by_pathology <- fread("/Users/Ewann/splicing_comparison/data/nygc/expression_by_pathology.csv")
annotation <- fread("/Users/Ewann/splicing_comparison/data/GRCh38/gencode.v42.annotation.gtf.genes.tab")

# looking at gene expression in NYGC vs detection rate for cryptic events 
potential_new_selective = expression_by_pathology |>
  filter(fraction_not_path <= 0.005) |>
  filter(fraction_path >= 0.01)

mean_tpm = tpm_nygc |> 
  mutate(gene = gsub("\\..*", "", gene)) |> 
  mutate(mean_tpm = rowMeans(select(tpm_nygc, contains("CGND")))) |> 
  select(gene, mean_tpm) |> 
  group_by(gene) |> 
  summarise(mean_tpm = mean(mean_tpm)) |> 
  ungroup()


detection_tpm = annotation |> 
  rename(id = `Gene ID`, gene = `Gene Symbol`) |> 
  select(id, gene) |> 
  mutate(id = gsub("\\..*", "", id)) |> 
  unique() |> 
  right_join(mean_tpm, by = c("id" = "gene")) |> 
  select(-id) |> 
  group_by(gene) |> 
  summarise(mean_tpm = mean(mean_tpm)) |> 
  ungroup() |> 
  unique() |> 
  right_join(potential_new_selective)

detection_tpm |> 
  ggplot(aes(x = fraction_path, y = mean_tpm)) +
  geom_point() +
  geom_smooth() +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.001)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1)) +
  ggrepel::geom_text_repel(aes(label = gene), size = 3) +
  theme_minimal()