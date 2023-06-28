library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

# read in nygc data
cell.type <- fread("/Users/Ewann/splicing_comparison/data/nygc/All_1917_samples_Darmanis_dtangle_deconv.tsv")
tpm_nygc = fread("/Users/Ewann/splicing_comparison/data/nygc/rsem_tpm_nygc.csv", header = TRUE)
nygc = arrow::read_parquet("/Users/Ewann/splicing_comparison/data/nygc/selective_cryptic_psi_in_nygc.parquet")
expression_by_pathology <- fread("/Users/Ewann/splicing_comparison/data/nygc/expression_by_pathology.csv")
annotation <- fread("/Users/Ewann/splicing_comparison/data/GRCh38/gencode.v42.annotation.gtf.genes.tab")


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

kalrn = c("chr3:124701255-124702038", "chr3:124701598-124702038", "chr3:124700033-124700975", "chr3:124700033-124701093")

nygc_psi = nygc |> 
  filter(tdp_path == "path") |>
  filter(paste_into_igv_junction %in% kalrn) |> 
  mutate(paste_into_igv_junction = gsub("[:-]", "_", paste_into_igv_junction)) |> 
  select(paste_into_igv_junction, sample, psi, disease) |> 
  pivot_wider(id_cols = c("sample", "disease"), names_from = "paste_into_igv_junction", values_from = "psi")
  

# add rbp and tdp columns
selective_gene <-data.table(
  id = c("ENSG00000104852", "ENSG00000112081", "ENSG00000204574", "ENSG00000160145", "ENSG00000120948", "ENSG00000116560", "ENSG00000037474"),
  gene = c("SNRNP70", "SRSF3", "ABCF1", "KALRN", "TARDBP", "SFPQ", "NSUN2")
)

tmp = nygc_by_gene |> 
  select(selective_gene$id, sample) |> 
  mutate(across(!sample, ~scale(.), .names = "{.col}"))
colnames(tmp) <- c(selective_gene$gene, "sample")

metatable = nygc_psi |> 
  left_join(tmp) |> 
  left_join(cell_type)  |> 
  tibble::column_to_rownames(var = "sample") |> 
  filter(grepl("Cortex", tissue_clean, ignore.case = TRUE))


# perform regression
plot(metatable |> select(-disease, -tissue_clean))

# KALRN MR modelling and model comparison
lm(ABCF1 ~ TARDBP + neurons + astrocytes + endothelial + microglia + disease, data = metatable) |> summary()

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
  select(-disease, -tissue_clean)) |> 
  tibble::column_to_rownames("rowname") |> 
  as.matrix()

cormat = cor(metatable |> select(-disease, -tissue_clean))

corrplot::corrplot(cormat, 
                   method = "color", 
                   type = "lower", 
                   tl.col = "black", 
                   tl.srt = 45, 
                   tl.cex = 0.8,
                   p.mat = pmat,
                   sig.level = 0.05,
                   insig = "blank")

# SRSF3 cell-type examination
metatable |> 
  tibble::rownames_to_column("sample") |> 
  pivot_longer(cols = c("neurons", "oligodendrocytes", "astrocytes", "microglia", "endothelial"), names_to = "cell_type", values_to = "deconv")  |>
  left_join(cell.type |> select(sample, tissue_clean)) |> 
  ggplot(aes(x = TARDBP, y = ABCF1, colour = disease)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  # ggpubr::stat_cor(size = 3) +
  # scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  # scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  # facet_grid(disease~tissue_clean) +
  theme_minimal()


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


# =====
metatable |> 
  ggplot(aes(x = KALRN, y = chr3_124701255_124702038)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor() +
  theme_minimal()
