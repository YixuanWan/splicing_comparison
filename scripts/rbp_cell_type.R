library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

# read in nygc data
cell.type <- fread("/Users/Ewann/splicing_comparison/data/nygc/All_1917_samples_Darmanis_dtangle_deconv.tsv")
tpm_nygc = fread("/Users/Ewann/splicing_comparison/data/nygc/rsem_tpm_nygc.csv", header = TRUE)
nygc = arrow::read_parquet("/Users/Ewann/splicing_comparison/data/nygc/selective_cryptic_psi_in_nygc.parquet")


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
  # filter(tdp_path == "path") |>
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
  left_join(cell_type) |> 
  column_to_rownames(var = "sample") 


# perform regression
plot(metatable |> select(-disease))

# KALRN MR modelling and model comparison
lm(chr3_124701255_124702038 ~ TARDBP + neurons + astrocytes + disease + tissue_clean, data = metatable) |> summary()

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
  select(-disease)) |> 
  column_to_rownames("rowname") |> 
  as.matrix()

cormat = cor(metatable |> select(-disease))

corrplot::corrplot(cormat, 
                   method = "color", 
                   type = "lower", 
                   tl.col = "black", 
                   tl.srt = 45, 
                   p.mat = pmat,
                   sig.level = 0.05,
                   insig = "blank")

# SRSF3 cell-type examination
metatable |> 
  rownames_to_column("sample") |> 
  pivot_longer(cols = c("neurons", "oligodendrocytes", "astrocytes", "microglia", "endothelial"), names_to = "cell_type", values_to = "deconv")  |>
  left_join(cell.type |> select(sample, tissue_clean)) |> 
  ggplot(aes(x = TARDBP, y = SNRNP70, colour = disease)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  facet_wrap(~tissue_clean, nrow = 2) 


