library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(stringr)
library(GenomicRanges)
library(rtracklayer)
library(pheatmap)
library(cluster)
library(dbscan)
library(ggcorrplot)


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
  filter(disease_tissue == TRUE) |> 
  select(paste_into_igv_junction, psi, tissue_clean, disease) 

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
  select(gene, junction, psi) |> 
  pivot_wider(names_from = c(gene, junction), 
              values_from = psi,
              values_fn = list) |> 
  purrr::map(as.data.table)

psi_wide = do.call(cbind, psi_tdp_wide) |> setnames(names(psi_tdp_wide))

# make psi_long for clustering
psi_long = t(psi_wide) |> as.data.frame(row.names = colnames(psi_wide)) 

# kmeans clustering - decide the number of clusters
# NOTE: substituted dataset with kmeans_higher_psi - the cluster with higher PSI in the original psi_long data
  
k <- 2:10

# silhouette scores
avg_sil = purrr::map(k, function(k){
                        km = kmeans(psi_long, centers = k, nstart = 25)
                        ss = silhouette(km$cluster, dist(psi_long))
                        mean(ss[, 3])
                        })
plot(k, type='b', avg_sil, 
     xlab='Number of clusters', 
     ylab='Average Silhouette Scores', 
     frame=FALSE) #k = 2


# clustering - actual clustering
kmeans = kmeans(psi_long, 2)

kmeans_psi = data.frame(mean = rowMeans(psi_long), cluster = kmeans$cluster)
kmeans_psi$gene = rownames(psi_long)


# plot the results with ggplot2
pca = prcomp(psi_long, center = TRUE)
pc_psi = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])

kmeans_psi_plot = cbind(kmeans_psi, pc_psi)
kmeans_psi_plot$cluster <- factor(kmeans_psi_plot$cluster)
kmeans_psi_plot |> 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(
    # fill = mean,
    # shape = cluster,
    colour = cluster), 
    # colour = "grey21",
    size = 2.3) +
  labs(title = "K-Means clustering results with average PSI",
       x = "PC1",
       y = 'PC2',
       colour =
       #   "average PSI",
       # shape = 
         "Cluster") +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 2)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 2)) +
  scale_fill_gradient(
    high = "turquoise",
    low = "salmon",
    space = "Lab",
    guide = "colourbar",
    aesthetics = "fill") +
  scale_shape_manual(values = c(21:25)) +
  theme_minimal()

# distribution of average PSI in K-means clusters
kmeans_psi_plot |> 
  ggplot(aes(x = cluster, y = rowMeans.psi_long.)) +
  geom_boxplot() +
  geom_jitter(aes(colour= cluster), size = 1) +
  labs(
    x = "K-means cluster",
    y = "Average PSI",
    colour = "Cluster"
  ) +
  theme_minimal()

# kmeans clustering with kmeans_higher_psi - the cluster with higher PSI in the original psi_long data
kmeans_higher_psi = psi_long |> 
  mutate(cluster = dbscan_psi_plot$cluster) |> 
  filter(cluster == 0) |> 
  select(-cluster)

k <- 2:10

# silhouette scores
avg_sil = purrr::map(k, function(k){
  km = kmeans(kmeans_higher_psi, centers = k, nstart = 25)
  ss = silhouette(km$cluster, dist(kmeans_higher_psi))
  mean(ss[, 3])
})
plot(k, type='b', avg_sil, 
     xlab='Number of clusters', 
     ylab='Average Silhouette Scores', 
     frame=FALSE) #k = 2


# clustering - actual clustering
kmeans = kmeans(kmeans_higher_psi, 2)

kmeans_psi = data.frame(mean = rowMeans(kmeans_higher_psi), cluster = kmeans$cluster)
kmeans_psi$gene = rownames(kmeans_higher_psi)


# plot the results with ggplot2
pca = prcomp(kmeans_higher_psi, center = TRUE)
pc_psi = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])

kmeans_psi_plot = cbind(kmeans_psi, pc_psi)
kmeans_psi_plot$cluster <- factor(kmeans_psi_plot$cluster)
kmeans_psi_plot |> 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = cluster), size = 2.3) +
  labs(title = "K-Means clustering results in cluster 1",
       x = "PC1",
       y = 'PC2',
       colour = "Cluster") +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 2)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 2)) +
  # scale_fill_gradient(
  #   high = "turquoise",
  #   low = "salmon",
  #   space = "Lab",
  #   guide = "colourbar",
  #   aesthetics = "fill") +
  # scale_shape_manual(values = c(21:25)) +
  theme_minimal()

# clustering with DBSCAN
dbscan_psi = dbscan(psi_long, eps = 1, minPts = 5)
dbscan_psi_plot = data.frame(rowMeans(psi_long), cluster = dbscan_psi$cluster)
dbscan_psi_plot$cluster = factor (dbscan_psi_plot$cluster)
dbscan_psi_plot$gene = rownames(psi_long)
dbscan_psi_plot = cbind(dbscan_psi_plot, pc_psi)

dbscan_psi_plot |> 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = rowMeans.psi_long., shape = cluster), colour = 'grey21', size = 2.3) +
  labs(title = "DBSCAN clustering results with average PSI",
       x = "PC1",
       y = 'PC2',
       colour = "average PSI",
       shape = "Cluster") +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  scale_fill_gradient(
    high = "turquoise",
    low = "salmon",
    space = "Lab",
    guide = "colourbar",
    aesthetics = "fill") +
  scale_shape_manual(values = c(21:25)) +
  theme_minimal()

# distribution of average PSI in dbscan clusters
dbscan_psi_plot |> 
  ggplot(aes(x = cluster, y = rowMeans.psi_long.)) +
  geom_boxplot() +
  geom_jitter(aes(colour= cluster), size = 1) +
  labs(
    x = "DBSCAN cluster",
    y = "Average PSI",
    colour = "Cluster"
  ) +
  theme_minimal()

# DBSCAN in cluster 0
higher_psi = psi_long |> 
  mutate(cluster = dbscan_psi_plot$cluster) |> 
  filter(cluster == 0) |> 
  select(-cluster)

dbscan_higher_psi = dbscan(higher_psi, eps = .3, minPts = 3)
dbscan_higher_psi_plot = data.frame(meanPSI = rowMeans(higher_psi), cluster = dbscan_higher_psi$cluster)
dbscan_higher_psi_plot$cluster = factor (dbscan_higher_psi_plot$cluster)
dbscan_higher_psi_plot$gene = rownames(higher_psi)

pca_higher = prcomp(higher_psi, center = TRUE)
pc_higher_psi = data.frame(PC1 = pca_higher$x[,1], PC2 = pca_higher$x[,2])

dbscan_higher_psi_plot = cbind(dbscan_higher_psi_plot, pc_higher_psi)

dbscan_higher_psi_plot |> 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = meanPSI, shape = cluster), colour = 'grey21', size = 2.3) +
  labs(title = "DBSCAN clustering results with average PSI",
       x = "PC1",
       y = 'PC2',
       colour = "average PSI",
       shape = "Cluster") +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  scale_fill_gradient(
    high = "turquoise",
    low = "salmon",
    space = "Lab",
    guide = "colourbar",
    aesthetics = "fill") +
  scale_shape_manual(values = c(21:25)) +
  theme_minimal()

# use pheatmap to plot the results
psi_sub = psi_long |> 
  mutate(cluster = dbscan_psi_plot$cluster) |> 
  filter(cluster == 0) |> 
  select(-cluster)

gene_cluster = cbind(kMeans = kmeans_psi$cluster, DBSCAN = dbscan_psi$cluster) |> 
  as.data.table() |> 
  mutate(kMeans = factor(kMeans), .keep = "unused") |> 
  mutate(DBSCAN = factor(DBSCAN), .keep = "unused")
rownames(gene_cluster) = cluster_psi$gene

kmeans_gene = as.data.table(kmeans_psi_plot$cluster) |> 
  mutate(cluster = factor(V1), .keep = "unused")  
rownames(kmeans_gene) = kmeans_psi_plot$gene

pheatmap(kmeans_higher_psi, 
         show_rownames = FALSE, 
         show_colnames = FALSE, 
         # cluster_cols = FALSE,
         annotation_row = kmeans_gene,
         fontsize = 5,
         main = 'junctions in K-means cluster 1')

hclust_kmeans = hclust(dist(kmeans_higher_psi), method = "complete")
stats::as.dendrogram(hclust_kmeans) |> 
  plot(horiz = TRUE)





# perform correlation
cor_matrix = psi_wide |> cor() |> round(3)
p_cor_matrix = psi_wide |> cor_pmat() 
padj_cor_matrix = p.adjust(p_cor_matrix, method = "BY", n = length(p_cor_matrix))  |> 
  matrix(ncol = dim(p_cor_matrix)[1], dimnames = dimnames(cor_matrix))


# plot correlation matrix with p-values
ggcorrplot(cor_matrix,
           # type = "upper",
           hc.order = TRUE,
           ggtheme = theme_minimal(),
           outline.color = 'white',
           p.mat = padj_cor_matrix,
           sig.level = 0.05,
           insig = "blank") +
  theme(axis.text.x = element_text(size = 1),
        axis.text.y = element_text(size = 1))
  

# subsetting cortex data and plot correlation matrix
psi_tdp_cortex = psi_gene_tdp|> 
  filter(grepl("Cortex", tissue_clean)) |> 
  select(gene, junction, psi) |> 
  pivot_wider(names_from = c(gene, junction), 
              values_from = psi,
              values_fn = list) |> 
  purrr::map(as.data.table)

sd_psi = sapply(psi_cortex, sd) 
psi_cortex = do.call(cbind, psi_tdp_cortex) |> 
  setnames(names(psi_tdp_cortex)) |> 
  select(which(sd_psi != 0)) #ENAH_chr1:225535563-225537167

cor_cortex = psi_cortex |> cor() |> round(3)
p_cor_cortex = psi_cortex |> cor_pmat() 
padj_cor_cortex = p.adjust(p_cor_cortex, method = "BY", n = length(p_cor_cortex)) |> 
  matrix(ncol = dim(p_cor_cortex)[1], dimnames = dimnames(cor_cortex))

ggcorrplot(cor_cortex,
           type = "upper",
           hc.order = TRUE,
           ggtheme = theme_minimal(),
           outline.color = 'white',
           p.mat = padj_cor_cortex,
           sig.level = 0.05,
           insig = "label") +
  theme(axis.text.x = element_text(size = 1),
        axis.text.y = element_text(size = 1))


# Convert the correlation matrix to a long format data frame
cor_df = as.data.frame(as.table(cor_cortex))
padj_df = as.data.frame(as.table(padj_cor_cortex))

cortex_cor = cor_df |> 
  left_join(padj_df, by = c("Var1" = "Var1", "Var2" = "Var2")) |> 
  mutate(estimate = Freq.x, .keep = "unused") |> 
  mutate(padj = Freq.y, .keep = "unused") |> 
  filter(Var1 != Var2) 

#get rid of the duplicates
df_sorted <- t(apply(cortex_cor, 1, function(x) sort(x)))
cor_unique = cortex_cor[!duplicated(df_sorted), ]

#get rid of the donor/accepter pairs
cor_gene1 = data.frame(str_split_fixed(cor_unique$Var1, "[_]", 2)) 
names(cor_gene1) <- c("gene1", "junction1")
cor_gene2 = data.frame(str_split_fixed(cor_unique$Var2, "[_]", 2)) 
names(cor_gene2) <- c("gene2", "junction2")

cor_gene = cbind(cor_unique, cor_gene1, cor_gene2) |>  
  select(-Var1, -Var2) |> 
  filter(gene1 != gene2) |> 
  arrange(desc(estimate), padj)

# plot the correlation
sig_cor_cortex = cor_gene |> 
  filter(padj < 0.05) 

# sig_cor_cortex |> 
#   ggplot(aes(x = gene1, y = gene2, fill = estimate)) +
#   geom_tile() + 
#   theme_minimal() +
#   theme(axis.text.x = element_text(size = 1),
#         axis.text.y = element_text(size = 1)) +
#   scale_fill_gradient(low = "white", high = "blue")



# write out the results
write.csv(cor_df, "/Users/Ewann/splicing_comparison/results/nygc/cryptic_tdp_corr.csv", row.names = FALSE)
write.csv(cor_gene, "/Users/Ewann/splicing_comparison/results/nygc/cryptic_tdp_corr_by_gene.csv", row.names = FALSE)
