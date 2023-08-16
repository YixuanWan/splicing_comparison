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

invitro_psi = splicingTable |> 
  filter(paste_into_igv_junction %in% kalrn) |> 
  select(gene_name,paste_into_igv_junction,junc_cat,baseline_PSI,contrast_PSI,comparison,strand, junc_cat) |>  
  mutate(is_cryptic = contrast_PSI > 0.1 & baseline_PSI < 0.05)  |> 
  filter(paste_into_igv_junction == "chr3:124701255-124702038") |> 
  group_by(comparison, paste_into_igv_junction) |> 
  summarise(contrast_PSI = mean(contrast_PSI), baseline_PSI = mean(baseline_PSI), is_cryptic) |> 
  unique() |> 
  right_join((experiments |> select(comparison, experiment, cell.type))) |> 
  filter(!is.na(paste_into_igv_junction)) 

invitro_psi  |> 
  separate(comparison, c("ctrl", "kd")) |> 
  # mutate(kd = forcats::fct_reorder(kd, contrast_PSI)) |> 
  # mutate(ctrl = forcats::fct_reorder(ctrl, baseline_PSI)) |>
  mutate(fc_PSI = contrast_PSI/baseline_PSI) |> 
  mutate(ctrl = forcats::fct_reorder(ctrl, fc_PSI)) |> 
  ggplot(aes(x = log2(fc_PSI), y = ctrl, fill = cell.type)) +
  geom_col() +
  geom_text(aes(label = experiment), size = 3, nudge_x = -1.5) +
  theme_minimal() + 
  # geom_vline(xintercept = 0.05, linetype = 2, colour = "turquoise") +
  theme(axis.text.y = element_blank()) +
  labs(
    title = "KALRN_chr3:124701255-124702038",
    subtitle = "PSI in in vitro KD studies",
    x = "log2(ContrastPSI/BaselinePSI)",
    y = "Experiment",
    fill = "Cell Type"
  ) 

cryptic_experiments = invitro_psi |> 
  pull(comparison)

# RBP cell-type differences check
filepath = list.files("/Users/Ewann/splicing_comparison/data/deseq2", full.names = TRUE)
suffix = ".csv"
filenames = list.files("/Users/Ewann/splicing_comparison/data/deseq2")
experiment_names = gsub(suffix,"",filenames)
estimate_files = purrr::map(filepath,my_clean_reader)
estimate_files = purrr::map2(estimate_files, experiment_names, ~cbind(.x, deseq2_table_name = .y))

RBP = c("SRSF3", "SNRNP70", "ABCF1", "KALRN", "NSUN2", "RBFOX2", "NUFIP2", "NPM1", "ELAVL3")
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

GENE = "NPM1"
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

# RBP correlation with number of cryptics

experiments = experiments |> 
  left_join(deseq2) 

experiments = experiments |> 
  mutate(kalrn_cryptic = comparison %in% cryptic_experiments)

GENE = "SNRNP70"
experiments |> 
  ggplot(aes(x = .data[[paste0("log2fold_change_", GENE)]], y = n_cryptic_junctions, colour = kalrn_cryptic)) +
  geom_point() +
  # ggforce::geom_mark_ellipse(data = (experiments |> filter(kalrn_cryptic == TRUE))) +
  geom_smooth(data = (experiments |> filter(kalrn_cryptic == TRUE)), method = "lm", size = 0.3) +
  ggpubr::stat_cor(data = (experiments |> filter(kalrn_cryptic == TRUE)), size = 3, label.x.npc = "middle") +
  # geom_vline(xintercept = 0, colour = 'grey', linetype = "dashed") +
  theme_minimal() +
  # ggrepel::geom_text_repel(data = (experiments |> filter(grepl("CHX", experiment))), aes(label = experiment), size = 2, colour = 'black') +
  labs(
    x = paste0("log2FC/", GENE),
    y = "N of cryptic junctions",
    colour = "KALRN cryptic?"
  ) +
  scale_colour_manual(values = c('lightblue', "darkorange"))

# cell-type variation in RBP levels
GENE = "tdp"
experiments |> 
  ggplot(aes(x = kalrn_cryptic, y = .data[[paste0("log2fold_change_", GENE)]], fill = kalrn_cryptic)) +
  geom_boxplot() +
  geom_jitter() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  ggpubr::stat_compare_means(size = 3) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    x = "With KALRN cryptic?",
    # y = paste0("log2FC/", GENE)
    y = "log2FC/TARDBP"
  ) +
  scale_fill_manual(values = c('lightblue', "orange"))

# correlation with tdp levels
experiments |> 
  ggplot(aes(y = .data[[paste0("log2fold_change_", GENE)]], x = log2fold_change_tdp, colour = kalrn_cryptic)) +
  geom_point() +
  # ggforce::geom_mark_ellipse(data = (experiments |> filter(kalrn_cryptic == TRUE))) +
  geom_smooth(data = (experiments |> filter(kalrn_cryptic == TRUE)), method = "lm", size = 0.3) +
  ggpubr::stat_cor(data = (experiments |> filter(kalrn_cryptic == TRUE)), size = 3) +
  # geom_vline(xintercept = 0, colour = 'grey', linetype = "dashed") +
  theme_minimal() +
  # ggrepel::geom_text_repel(aes(label = experiment), size = 2) +
  labs(
    y = paste0("log2FC/", GENE),
    x = "log2FC/TARDBP",
    colour = "KALRN cryptic?"
  ) +
  scale_colour_manual(values = c('lightblue', "darkorange"))
