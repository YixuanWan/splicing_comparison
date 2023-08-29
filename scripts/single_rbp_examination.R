library(dplyr)
library(tidyr)
library(data.table)
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
experiments = read.csv("/Users/Ewann/splicing_comparison/samplesheet/tdp_experiments_updated-tdp_experiments_short.csv")

kalrn = c("chr3:124701255-124702038", "chr3:124701598-124702038", "chr3:124700033-124700975", "chr3:124700033-124701093")

invitro_psi = splicingTable |> 
  # filter(paste_into_igv_junction %in% kalrn) |>
  # select(gene_name,paste_into_igv_junction,junc_cat,baseline_PSI,contrast_PSI,comparison,strand, junc_cat) |>  
  mutate(is_cryptic = contrast_PSI > 0.1 & baseline_PSI < 0.05)  |>
  # filter(grepl("chr19:49101471-491021", paste_into_igv_junction)) |> 
  filter(paste_into_igv_junction == "chr19:49101471-49102114") |>
  group_by(comparison, paste_into_igv_junction) |> View()
  summarise(contrast_PSI = mean(contrast_PSI), baseline_PSI = mean(baseline_PSI), is_cryptic) |> 
  unique() |> 
  right_join((experiments |> select(comparison, experiment, cell.type))) |> 
  filter(!is.na(paste_into_igv_junction)) 

invitro_psi  |> 
  # filter(cell.type == "SH-SY5Y") |> 
  # separate(comparison, c("ctrl", "kd")) |>
  # mutate(kd = gsub("tdp43kd","", kd)) |> 
  # mutate(kd = forcats::fct_reorder(kd, contrast_PSI)) |> 
  pivot_longer(cols = contains("PSI"), names_to = "condition", values_to = "psi") |>
  ggplot(aes(x = psi, y = comparison)) +
  geom_col(aes(fill = condition), position = 'dodge') +
  theme_minimal() + 
  # geom_vline(xintercept = 0.05, linetype = 2, colour = "#00BFC4") +
  # geom_vline(xintercept = 0.10, linetype = 2, colour = "#F8766D") +
  labs(
    title = "SNRNP70 AS",
    # subtitle = "CHX SH-SY5Y",
    x = "PSI",
    y = "Experiment",
    fill = "Condition"
  ) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"), labels = c("Control", "TDP-43 KD"))

cryptic_experiments = invitro_psi |> 
  pull(comparison)

# RBP cell-type differences check
filepath = list.files("/Users/Ewann/splicing_comparison/data/deseq2", full.names = TRUE)
suffix = ".csv"
filenames = list.files("/Users/Ewann/splicing_comparison/data/deseq2")
experiment_names = gsub(suffix,"",filenames)
estimate_files = purrr::map(filepath,my_clean_reader)
estimate_files = purrr::map2(estimate_files, experiment_names, ~cbind(.x, deseq2_table_name = .y))

RBP = c("SNRNP70", "KALRN", "ELAVL3", "SRSF3")
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

GENE = "SNRNP70"
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

GENE = "SRSF3"
experiments |> 
  ggplot(aes(x = .data[[paste0("log2fold_change_", GENE)]], y = n_cryptic_junctions)) +
  geom_smooth(method = 'lm', colour = "black", linewidth = 0.5) +
  geom_point(aes(colour = cell.type)) +
  # ggforce::geom_mark_ellipse(data = (experiments |> filter(kalrn_cryptic == TRUE))) +
  ggpubr::stat_cor(size = 3, label.x.npc = "left", label.y = 1050) +
  # geom_vline(xintercept = 0, colour = 'grey', linetype = "dashed") +
  theme_minimal() +
  # ggrepel::geom_text_repel(aes(label = experiment), size = 3, colour = 'black') +
  theme(legend.position = "none") +
  labs(
    x = paste0("log2FC/", GENE),
    y = "N of cryptic junctions"
  ) +
  ylim(0, 1200)

# cell-type variation in RBP levels across experiments
experiments |> 
  ggplot(aes(x = cell.type, y = .data[[paste0("log2fold_change_", GENE)]])) +
  geom_boxplot(colour = "lightgrey") +
  geom_jitter(aes(colour = cell.type)) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggpubr::stat_compare_means(method = "wilcox.test",
                             comparison = list(c(1,2), c(2,3), c(3,4), c(1,3),  c(2,4), c(1,4)),
                             label.y = c(0.4, 0.4, 0.4, 0.7, 0.5, 0.6), 
                             label = "p.format",
                             size = 3) +
  labs(
    x ="Cell Type",
    y = paste0("log2FC/", GENE))

# cell-type variation in RBP levels in experiments
experiments |> 
  mutate(comparison = glue::glue("{experiment}_{comparison}")) |> 
  mutate(log2fold_change_TARDBP = log2fold_change_tdp, .keep = 'unused') |> 
  pivot_longer(cols = contains("log2fold_change"), names_to = "RBP", names_prefix = "log2fold_change_", values_to = "log2FC") |>  
  select(-contains("padj")) |> 
  filter(RBP != "KALRN") |> 
  ggplot(aes(x = comparison, y = log2FC, fill = experiment)) +
  # geom_col() +
  # coord_flip() +
  # geom_boxplot(colour = 'black') +
  # geom_jitter() +
  # geom_hline(yintercept = 0, linetype = 'dashed') +
  # ggpubr::stat_compare_means(size = 3, label = "p.format") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.position = "none") +
  # labs(
  #   x = "KALRN SJ1?",
  #   y = "log2FC"
  # ) +
  # scale_colour_manual(values = c('lightblue', "orange")) +
  scale_fill_manual(values = c("#CC79A7", "#0072B2")) +
  facet_wrap(~RBP, nrow = 1)

# correlation with tdp levels
GENE = "SNRNP70"
experiments |> 
  ggplot(aes(y = .data[[paste0("log2fold_change_", GENE)]], x = log2fold_change_tdp, colour = kalrn_cryptic)) +
  ggforce::geom_mark_hull(aes(fill = kalrn_cryptic), colour = "white", alpha = 0.1) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  # geom_smooth(data = (experiments |> filter(kalrn_cryptic == TRUE)), method = "lm", size = 0.3) +
  ggpubr::stat_cor() +
  # geom_vline(xintercept = 0, colour = 'grey', linetype = "dashed") +
  theme_minimal() +
  # ggrepel::geom_text_repel(aes(label = experiment), size = 2) +
  labs(
    y = paste0("log2FC/", GENE),
    x = "log2FC/TARDBP",
    colour = "KALRN SJ1?",
    fill = "KALRN SJ1?"
  ) +
  scale_colour_manual(values = c('lightblue', "darkorange")) +
  scale_fill_manual(values = c('lightblue', "darkorange"))

# RBP correlation with KALRN cryptic PSI
experiments |> 
  select(-contains("padj")) |> 
  mutate(log2fold_change_TARDBP = log2fold_change_tdp, .keep = "unused") |> 
  right_join(invitro_psi) |>
  mutate(fc_psi = log2(contrast_PSI/baseline_PSI)) |> 
  pivot_longer(cols = contains("log2"), names_to = "rbp", values_to = "log2fc", names_prefix = "log2fold_change_") |> 
  ggplot(aes(x = log2fc, y = contrast_PSI)) +
  geom_point(aes(colour = cell.type), size = 2) +
  # geom_smooth(method = "lm") +
  # ggpubr::stat_cor() +
  # ggrepel::geom_text_repel(aes(label = experiment), size = 3, position = "stack", max.overlaps = 2) +
  theme_minimal() +
  facet_wrap(~ rbp, nrow = 1) +
  labs(
    x = "log2FC",
    y = "log2(ContrastPSI/baselinePSI)",
    colour = "Cell Type"
  )


# KALRN cryptic - RBP association in facs neurons
splicingTable |> 
  select(gene_name,paste_into_igv_junction,junc_cat,baseline_PSI,contrast_PSI,comparison,strand, junc_cat) |>  
  filter(paste_into_igv_junction == "chr3:124701255-124702038") |> 
  filter(comparison == "controlliufacsneurons-tdp43kdliufacsneurons") |> 
  group_by(paste_into_igv_junction) |> 
  summarise(contrast_PSI = mean(contrast_PSI), baseline_PSI = mean(baseline_PSI)) |> 
  unique() 
