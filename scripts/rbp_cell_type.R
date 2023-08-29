library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

# read in nygc data
cell.type <- fread("/Users/Ewann/splicing_comparison/data/nygc/All_1917_samples_Darmanis_dtangle_deconv.tsv")
tpm_nygc = fread("/Users/Ewann/splicing_comparison/data/nygc/rsem_tpm_nygc.csv", header = TRUE)
nygc = arrow::read_parquet("/Users/Ewann/splicing_comparison/data/nygc/all_psi_in_nygc.parquet")



nygc_by_gene = tpm_nygc |> 
  mutate(gene = gsub("\\..*", "", gene)) |> 
  distinct(gene, .keep_all = TRUE) |> 
  tibble::column_to_rownames(var = "gene") |> 
  t() |> as.data.table(keep.rownames = TRUE) |> 
  rename(sample = rn)

cell_type = cell.type |> 
  select(sample, cell, deconv, tissue_clean) |> 
  unique() |> 
  pivot_wider(id_cols = c("sample", "tissue_clean"), names_from = "cell", values_from = "deconv") 

junction = c("chr3:124701255-124702038", 
             "chr3:124701598-124702038", 
             "chr3:124700033-124700975", 
             "chr3:124700033-124701093", 
             "chr19:49101471-49102114",
             "chr8:79611214-79616822")

nygc_psi = nygc |> 
  filter(paste_into_igv_junction %in% junction) |> 
  mutate(paste_into_igv_junction = gsub("[:-]", "_", paste_into_igv_junction)) |> 
  select(paste_into_igv_junction, sample, psi, disease, tdp_path, onset) |> 
  pivot_wider(id_cols = c("sample", "disease", "tdp_path", "onset"), names_from = "paste_into_igv_junction", values_from = "psi")
  

nygc_psi |> 
  # mutate(SNRNP70 = ifelse(chr19_49101471_49102114 > 0, TRUE, FALSE)) |>
  # mutate(KALRN = ifelse(chr3_124701255_124702038 > 0, TRUE, FALSE)) |> 
  ggplot(aes(x = chr8_79611214_79616822, y = chr19_49101471_49102114)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.y = element_blank()) +
  labs(
    x = "STMN2/PSI",
    y = "SNRNP70 AS/PSI"
  )


  

# add rbp and tdp columns
selective_gene <-data.table(
  id = c("ENSG00000104852", "ENSG00000112081", "ENSG00000160145", "ENSG00000120948", "ENSG00000196361"),
  gene = c("SNRNP70", "SRSF3", "KALRN", "TARDBP", "ELAVL3")
)

tmp = nygc_by_gene |> 
  select(selective_gene$id, sample) |> 
  mutate(across(!sample, ~scale(.), .names = "{.col}"))
colnames(tmp) <- c(selective_gene$gene, "sample")

metatable = nygc_psi |> 
  left_join(tmp) |> 
  left_join(cell_type)  |> 
  tibble::column_to_rownames(var = "sample")




# perform regression
plot(sc_metatable |> select(-disease, -tissue_clean))

# KALRN MR modelling and model comparison
lm(chr3_124701255_124702038 ~ TARDBP + SNRNP70 + ELAVL3 + disease + tissue_clean + onset, data = metatable) |> summary()

# correlation by disease type
metatable |> 
  ggplot(aes(x = neurons, y = metatable[["chr3:124701255-124702038"]], colour = disease)) +
  geom_point() +
  geom_smooth(method = "lm") +
  # scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  # scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  theme_minimal() +
  labs(
    y = "KALRN_chr3:chr3:124701255-124702038"
  )

# correlation matrix
pmat = rstatix::cor_pmat(metatable |> 
  select(-disease, -tissue_clean,-tdp_path, -SRSF3, -onset)) |> 
  tibble::column_to_rownames("rowname") |> 
  as.matrix()

cormat = cor(metatable |> select(-disease, -tissue_clean,-tdp_path, -SRSF3, -m))

corrplot::corrplot(cormat, 
                   method = "color", 
                   type = "lower", 
                   tl.col = "black", 
                   tl.srt = 45, 
                   tl.cex = 0.8,
                   p.mat = pmat,
                   sig.level = 0.05,
                   insig = "blank")

# RBP cell-type, disease-type, tissue-type examination
plot = metatable |> 
  tibble::rownames_to_column("sample") |> 
  pivot_longer(cols = c("neurons", "oligodendrocytes", "astrocytes", "microglia", "endothelial"), names_to = "cell_type", values_to = "deconv")  |>
  left_join(cell.type |> select(sample, tissue_clean)) 
  # filter(tdp_path == "path") |> 
  # filter(!tissue_clean %in% c("Occipital_Cortex", "Sensory_Cortex")) 

plot |> 
  filter(chr3_124701255_124702038 > 0) |>
  # filter(TARDBP >= 0 & KALRN >= 0) |>
  # filter(!grepl("non", disease)) |>
  filter(grepl("Cortex", tissue_clean)) |> 
  filter(!tissue_clean %in% c("Occipital_Cortex", "Sensory_Cortex")) |> 
  ggplot(aes(x = deconv, y = chr3_124701255_124702038, colour = disease)) +
  geom_point(size = 1) +
  # ggforce::geom_mark_hull(alpha = 0.3) +
  # geom_smooth(method = 'lm') +
  # ggpubr::stat_cor(data = (plot |> filter(disease == "FTD-TDP")), size = 3) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  facet_grid(tissue_clean ~ cell_type) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    y = "KALRN SJ1/PSI",
    x = "Deconv",
    colour = "Disease",
    fill = "Disease") +
# scale_fill_manual(values = c("#F8766D", "#A3A500", "#00BA38", "#E76BF3")) +
scale_colour_manual(values = c("#F8766D", "#A3A500", "#00BA38", "#E76BF3"))

# KALRN cryptic vs KALRN TPM in all samples
plot |> 
  # filter(chr3_124701255_124702038 > 0) |>
  # filter(grepl("Cortex", tissue_clean, ignore.case = TRUE)) |>
  # filter(!tissue_clean %in% c("Occipital_Cortex", "Sensory_Cortex")) |>
  # filter(tdp_path == "path") |>
  ggplot(aes(x = chr19_49101471_49102114, y = SNRNP70)) +
  geom_point(aes(colour = disease)) +
  geom_smooth(method = 'lm', colour = "black", linewidth = 1) +
  ggpubr::stat_cor(method = 'spearman', size = 3) +
  labs(
    x = "SNRNP70 AS/PSI",
    y = "SNRNP70/TPM",
    colour = "") +
  # scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.01, base = 2)) +
  theme_minimal() +
  # scale_colour_manual(values = c("#F8766D", "#A3A500", "#00BA38", "#E76BF3")) +
  theme(legend.position = "bottom") 
  # facet_wrap(~tissue_clean)
  
  

# splicing isoform distribution by disease type
metatable |> 
  filter(grepl("Cortex", tissue_clean)) |>
  # filter(chr3_124701255_124702038 > 0) |>
  # filter(KALRN >= 0) |>
  # filter(!tissue_clean %in% c("Occipital_Cortex", "Sensory_Cortex")) |> 
  ggplot(aes(x = tdp_path, y = chr19_49101471_49102114, fill = tdp_path)) +
  geom_boxplot() +
  # geom_jitter(size = 0.5, colour = "grey") +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             # comparison = list(c(1,2), c(2,3), c(3,4)),
                             label = "p.format") +
                             # size = 2.8,
                             # label.y = c(-0.2, 0, 0)) +
  theme_minimal() +
  facet_wrap(~tissue_clean, nrow = 1) +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom") +
  labs(fill = "",
       x = "",
       y = "SNRNP70 AS/PSI") 
  # scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2))
  # ylim(-1,0.2) +
  # scale_fill_manual(values = c("#F8766D", "#A3A500", "#00BA38", "#E76BF3"))


# rbp ~ KALRN
metatable |> 
  # filter(grepl("ALS", disease)) |>
  filter(grepl("Cortex", tissue_clean, ignore.case = TRUE)) |>
  filter(!tissue_clean %in% c("Occipital_Cortex", "Sensory_Cortex")) |>
  # filter(tissue_clean %in% c("Frontal_Cortex", "Temporal_Cortex")) |>
  pivot_longer(cols = c("TARDBP", "SNRNP70", "ELAVL3"), names_to = "RBP", values_to = "TPM") |> 
  # filter(chr3_124701255_124702038 > 0) |>
  # filter(TPM >= 0) |>
  # filter(KALRN >= 0) |>
  ggplot(aes(x = disease, y = TPM, fill = disease)) +
  # geom_smooth(method = "lm", colour = 'black', linewidth = 0.8) +
  # geom_point() +
  # ggpubr::stat_cor(method = "spearman", size = 2.5, cor.coef.name = "rho", label.x.npc = "middle") +
  geom_boxplot() +
  ggpubr::stat_compare_means(method = "wilcox.test",
                             comparison = list(c(1,2), c(2,3), c(4,5), c(3,5)),
                             label = "p.signif",
                             size = 3,
                             label.y = c(3, 4.5, 3, 4.6)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_blank()) +
  facet_grid(tissue_clean ~ RBP) +
  # facet_wrap(~RBP, nrow = 1) +
  # scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  labs(x = "",
       y = "RBP/TPM",
       fill = "") 
  # scale_colour_manual(values = c("#F8766D", "#A3A500", "#00BA38", "#E76BF3"))
  # scale_color_manual(values = c("#A3A500", "#00BA38", "#E76BF3"))

# KALRN SJ1 containing samples vs others
metatable |> 
  # filter(!grepl("non", disease)) |>
  filter(grepl("Cortex", tissue_clean, ignore.case = TRUE)) |> 
  filter(!tissue_clean %in% c("Occipital_Cortex", "Sensory_Cortex")) |>
  pivot_longer(cols = c("TARDBP", "SNRNP70", "ELAVL3"), names_to = "RBP", values_to = "TPM") |> 
  mutate(cryptic = ifelse(chr3_124701255_124702038 > 0, TRUE, FALSE)) |> 
  ggplot(aes(x = TPM, y = KALRN, colour = chr19_49101471_49102114)) +
  geom_point(size = 1) +
  # geom_boxplot() +
  # scale_colour_manual(values = c('lightblue', "darkorange")) +
  # scale_fill_manual(values = c("#A3A500", "#00BA38", "#E76BF3")) +
  theme_minimal() +
  # ggpubr::stat_compare_means(method = "wilcox.test",
  #                            comparisons = list(c("ALS-TDP", "Control"), c("FTD-TDP", "Control"), c("ALS-TDP", "FTD-TDP")),
  #                            label = "p.signif",
  #                            size = 3,
  #                            # label.y = 5) +
  #                            label.y = c(3, 3.5, 5)) +
  facet_grid(tissue_clean ~ RBP) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank()) +
  labs(
    x = "RBP/TPM",
    y = "KALRN/TPM",
    colour = "SNRNP70 AS/PSI") +
  scale_colour_continuous(high = "#0072B2", low = "#CC79A7") +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2))



# =====SNRNP70 AS distribution by disease and tissue type
metatable |> 
  # filter(!grepl("non", disease)) |>
  filter(grepl("Cortex", tissue_clean)) |> 
  # filter(chr3_124701255_124702038 >= 0) |>
  filter(chr19_49101471_49102114 >= 0) |>
  # filter(SNRNP70 > 0) |>
  ggplot(aes(x = disease, y = chr19_49101471_49102114, fill = disease)) +
  # geom_point() +
  # geom_smooth(method = "lm") +
  geom_boxplot() +
  # ggpubr::stat_cor() +
  ggpubr::stat_compare_means(method = "wilcox.test",
                             comparison = list(c(1,2), c(2,3), c(4,5), c(3,5)),
                             label = "p.signif",
                             size = 3) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank()
  ) +
  facet_wrap(~tissue_clean, nrow = 1) +
  labs(
    y = "SNRNP70 AS/PSI",
    x = "Disease"
  ) 
  scale_color_manual(values = c("#A3A500", "#00BA38", "#E76BF3"))
  

# KALRN SJ1 containing samples, RBP level examination

  
