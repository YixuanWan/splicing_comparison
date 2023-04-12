library(arrow)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(stringr)
library(GenomicRanges)
library(rtracklayer)


# load in dataset
nygc = read_parquet("/Users/Ewann/splicing_comparison/data/nygc/selective_cryptic_psi_in_nygc.parquet")
gencode = as.data.table(import.bed("/Users/Ewann/splicing_comparison/data/GRCh38/gencode.v42.annotation.bed12"))

cluster_mount_point = "/Users/Ewann/cluster"
tx2gene = file.path(cluster_mount_point,"vyplab_reference_genomes/annotation/human/GRCh38/gencode.v42.tx2gene.csv")
tx2gene <- fread(tx2gene,header=FALSE)
colnames(tx2gene) = c("TXNAME", "GENEID")

# generate a seq2gene table
gencodes = gencode |> 
  left_join(tx2gene, by = c("name" = "TXNAME")) |> 
  left_join(gentab, by = c("GENEID" = "Gene.ID")) |> 
  select(seqnames, start, end, Gene.Symbol) |> 
  mutate(gene = Gene.Symbol, .keep = "unused")

gr_gencodes = GRanges(gencodes) |> unique()

# select tdp-pathology patients, get rid of other columns
nygc_tdp = nygc |>
  filter(tdp_path == "path") |> 
  select(paste_into_igv_junction, psi) 

# create a GRange list for psi_tdp
gene_tdp = data.frame(str_split_fixed(nygc_tdp$paste_into_igv_junction, "[:-]", 3)) 
names(gene_tdp) <- c("seqnames", "start", "end")
psi_tdp = cbind(nygc_tdp, gene_tdp) |> select(-paste_into_igv_junction)
gr_psi_tdp = GRanges(psi_tdp)

# get gene names for psi_tdp
overlaps_gene_tdp = findOverlaps(gr_psi_tdp, gr_gencodes)
uniqueHits = overlaps_gene_tdp |> 
  as.data.frame() |> 
  distinct(queryHits, .keep_all = TRUE)
overlaps_tdp = gr_psi_tdp[uniqueHits$queryHits]
mcols(overlaps_tdp)$gene = mcols(gr_gencodes[uniqueHits$subjectHits])$gene

psi_gene_tdp = overlaps_tdp |> 
  as.data.table() |> 
  mutate(junction = paste0(seqnames, ":", start, "-", end), .keep = "unused") |> 
  select(-width, -strand)

# make psi_wide for correlation 
psi_tdp_wide = psi_gene_tdp|> 
  pivot_wider(names_from = c(gene, junction), 
              values_from = psi,
              values_fn = list) |> 
  purrr::map(as.data.table)

psi_wide = do.call(cbind, psi_tdp_wide) |> setnames(names(psi_tdp_wide))

# perform correlation
cor_matrix = psi_wide |> cor() 
diag(cor_matrix) = 0

# Convert the correlation matrix to a long format data frame
cor_df = as.data.frame(as.table(cor_matrix))|> 
  arrange(desc(Freq)) 

cor_gene1 = data.frame(str_split_fixed(cor_df$Var1, "[_]", 2)) 
names(cor_gene1) <- c("gene1", "junction1")
cor_gene2 = data.frame(str_split_fixed(cor_df$Var2, "[_]", 2)) 
names(cor_gene2) <- c("gene2", "junction2")

cor_gene = cbind(cor_df, cor_gene1, cor_gene2) |> 
  select(-Var1, -Var2) |> 
  filter(gene1 != gene2) |> 
  arrange(desc(Freq))

# write out the results
write.csv(cor_df, "/Users/Ewann/splicing_comparison/results/nygc/cryptic_tdp_corr.csv", row.names = FALSE)
write.csv(cor_gene, "/Users/Ewann/splicing_comparison/results/nygc/cryptic_tdp_corr_by_gene.csv", row.names = FALSE)


# use pheatmap to plot the results
library(pheatmap)
cor_gene = read.csv("/Users/Ewann/splicing_comparison/results/nygc/cryptic_tdp_corr_by_gene.csv", header = TRUE)
pheatmap(cor_gene)

# Create the heatmap using ggplot2
top_cor |>   
ggplot(aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "lightgrey", high = "purple") +
  labs(x = "", y = "", fill = "Correlation") +
  theme_minimal() +
  # theme(axis.text.x = element_blank(), axis.text.y = element_blank())
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
