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

kalrn = c("chr3:124701255-124702038", "chr3:124701598-124702038", "chr3:124700033-124700975", "chr3:124700033-124701093", "chr19:49101471-49102114")

nygc_psi = nygc |> 
  filter(paste_into_igv_junction %in% kalrn) |> 
  mutate(paste_into_igv_junction = gsub("[:-]", "_", paste_into_igv_junction)) |> 
  select(paste_into_igv_junction, sample, psi, disease, tdp_path) |> 
  pivot_wider(id_cols = c("sample", "disease", "tdp_path"), names_from = "paste_into_igv_junction", values_from = "psi")
  

nygc_psi |> 
  mutate(SNRNP70 = ifelse(chr19_49101471_49102114 > 0, TRUE, FALSE)) |>
  mutate(KALRN = ifelse(chr3_124701255_124702038 > 0, TRUE, FALSE)) |> 
  ggplot(aes(x = disease, y = KALRN)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.y = element_blank()) +
  labs(
    x = "Disease",
    y = ""
  )


  

# add rbp and tdp columns
selective_gene <-data.table(
  id = c("ENSG00000104852", "ENSG00000112081", "ENSG00000204574", "ENSG00000160145", "ENSG00000120948", "ENSG00000116560", "ENSG00000037474", "ENSG00000196361"),
  gene = c("SNRNP70", "SRSF3", "ABCF1", "KALRN", "TARDBP", "SFPQ", "NSUN2", "ELAVL3")
)

tmp = nygc_by_gene |> 
  select(selective_gene$id, sample) |> 
  mutate(across(!sample, ~scale(.), .names = "{.col}"))
colnames(tmp) <- c(selective_gene$gene, "sample")

metatable = nygc_psi |> 
  left_join(tmp) |> 
  left_join(cell_type)  |> 
  tibble::column_to_rownames(var = "sample")

cortex_metatable = nygc_psi |> 
  left_join(tmp) |> 
  left_join(cell_type)  |> 
  tibble::column_to_rownames(var = "sample") |> 
  filter(grepl("Cortex", tissue_clean, ignore.case = TRUE))  



# perform regression
plot(sc_metatable |> select(-disease, -tissue_clean))

# KALRN MR modelling and model comparison
lm(chr3_124701255_124702038 ~ TARDBP + SNRNP70 + ELAVL3 + disease + tissue_clean, data = cortex_metatable) |> summary()

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
pmat = rstatix::cor_pmat(cortex_metatable |> 
  select(-disease, -tissue_clean, -SFPQ, -ABCF1, -NSUN2, -tdp_path)) |> 
  tibble::column_to_rownames("rowname") |> 
  as.matrix()

cormat = cor(cortex_metatable |> select(-disease, -tissue_clean, -SFPQ, -ABCF1, -NSUN2, -tdp_path))

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
plot = cortex_metatable |> 
  tibble::rownames_to_column("sample") |> 
  pivot_longer(cols = c("neurons", "oligodendrocytes", "astrocytes", "microglia", "endothelial"), names_to = "cell_type", values_to = "deconv")  |>
  left_join(cell.type |> select(sample, tissue_clean)) |> 
  # filter(tdp_path == "path") |> 
  filter(!tissue_clean %in% c("Occipital_Cortex", "Sensory_Cortex")) 

plot |> 
  # filter(chr3_124701255_124702038 > 0) |>
  filter(TARDBP >= 0 & KALRN >= 0) |> 
  filter(!grepl("non", disease)) |>
  ggplot(aes(x = TARDBP, y = KALRN, colour = disease)) +
  geom_point(size = 0.8) +
  ggforce::geom_mark_ellipse() +
  # geom_smooth(method = 'lm') +
  # ggpubr::stat_cor(data = (plot |> filter(disease == "FTD-TDP")), size = 3) +
  # scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  # scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  facet_wrap(~tissue_clean) +
  theme_minimal() +
  labs(
    colour = "Disease Type") +
  scale_color_manual(values = c("#A3A500", "#00BA38", "#E76BF3"))

# disease type
cortex_metatable |> 
  ggplot(aes(x = disease, y = KALRN, fill = disease)) +
  geom_boxplot() +
  ggpubr::stat_anova_test(size = 2.5, label.x.npc = "middle") +
  theme_minimal() +
  facet_wrap(~tissue_clean, nrow = 1) +
  theme(axis.text.x = element_blank()) +
  labs(fill = "Disease Type")

# tdp vs rbp ~ KALRN
cortex_metatable |> 
  filter(!grepl("non", disease)) |>
  filter(grepl("Cortex", tissue_clean, ignore.case = TRUE)) |>
  filter(!tissue_clean %in% c("Occipital_Cortex", "Sensory_Cortex")) |> 
  pivot_longer(cols = c("TARDBP", "SNRNP70"), names_to = "RBP", values_to = "TPM") |> 
  filter(chr19_49101471_49102114 >= 0) |>
  filter(TPM >= 0) |> 
  # filter(KALRN >= 0) |> 
  ggplot(aes(x = TPM, y = chr19_49101471_49102114, colour = disease)) +
  geom_point(size = 0.8) +
  geom_smooth(method = "lm") +
  # ggforce::geom_mark_ellipse() +
  # ggpubr::stat_cor(size = 3) +
  theme_minimal() +
  facet_grid(tissue_clean ~ RBP) +
  labs(x = "RBP") +
  scale_color_manual(values = c("#A3A500", "#00BA38", "#E76BF3"))
  # scale_color_brewer(palette = "Set1")


# =====
cortex_metatable |> 
  filter(!grepl("non", disease)) |>
  # filter(chr3_124701255_124702038 >= 0) |>
  filter(chr19_49101471_49102114 >= 0) |>
  filter(SNRNP70 > 0) |>
  ggplot(aes(x = SNRNP70, y = chr19_49101471_49102114, colour = disease)) +
  geom_point() +
  geom_smooth(method = "lm") +
  # ggpubr::stat_cor() +
  theme_minimal() +
  scale_color_manual(values = c("#A3A500", "#00BA38", "#E76BF3"))
