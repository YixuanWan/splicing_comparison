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

# make psi_long for clustering
psi_long = t(psi_wide) |> as.data.frame(row.names = colnames(psi_wide)) 

# clustering - decide the number of clusters
  
k <- 2:20

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

# elbow - wss
wss = purrr::map(k, function(k){
                    kmeans(psi_long, centers = k, nstart = 25)$tot.withinss
                  })
plot(k, type='b', wss, #plot wss
     xlab='Number of clusters', 
     ylab='WSS', 
     frame=FALSE) #k = 2


# clustering - actual clustering
kmeans_psi = kmeans(psi_long, 3)

cluster_psi = data.frame(rowMeans(psi_long), cluster = kmeans_psi$cluster)
cluster_psi$gene = rownames(psi_long)


# plot the results with ggplot2
pca = prcomp(psi_long)
pc_psi = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])

cluster_psi_plot = cbind(cluster_psi, pc_psi)
cluster_psi_plot$cluster <- factor(cluster_psi_plot$cluster)
cluster_psi_plot |> 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = cluster)) +
  labs(title = "K-means Clustering Results",
       x = "PC1",
       y = 'PC2',
       colour = "Cluster") +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  theme_minimal()

# clustering with DBSCAN
dbscan_psi = dbscan(psi_long, eps = 1, minPts = 5)
dbscan_psi_plot = data.frame(mean_psi, cluster = dbscan_psi$cluster)
dbscan_psi_plot$cluster = factor (dbscan_psi_plot$cluster)
dbscan_psi_plot$gene = rownames(psi_long)
dbscan_psi_plot = cbind(dbscan_psi_plot, pc_psi)

dbscan_psi_plot |> 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = cluster)) +
  labs(title = "DBSCAN Clustering Results",
       x = "PC1",
       y = 'PC2',
       colour = "Cluster") +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 10)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 10)) +
  theme_minimal()


# use pheatmap to plot the results
psi_sub = psi_long |> 
  filter(rowMeans(psi_long) > 0.05)

gene_cluster = cbind(kMeans = cluster_psi$cluster, DBSCAN = dbscan_psi$cluster) |> 
  as.data.table() |> 
  mutate(kMeans = factor(kMeans), .keep = "unused") |> 
  mutate(DBSCAN = factor(DBSCAN), .keep = "unused")
rownames(gene_cluster) = cluster_psi$gene

pheatmap(psi_sub, 
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         cluster_cols = FALSE,
         annotation_row = gene_cluster,
         cutree_rows = 2,
         fontsize = 5,
         main = 'junctions with average psi > 0.07')




# perform correlation
cor_matrix = psi_wide |> cor() 

# Convert the correlation matrix to a long format data frame
cor_df = as.data.frame(as.table(cor_matrix))|> 
  arrange(desc(Freq)) |> 
  filter(Var1 != Var2)

#get rid of the duplicates
df_sorted <- t(apply(cor_df, 1, function(x) sort(x)))
cor_unique = cor_df[!duplicated(df_sorted), ]

cor_gene1 = data.frame(str_split_fixed(cor_unique$Var1, "[_]", 2)) 
names(cor_gene1) <- c("gene1", "junction1")
cor_gene2 = data.frame(str_split_fixed(cor_unique$Var2, "[_]", 2)) 
names(cor_gene2) <- c("gene2", "junction2")

cor_gene = cbind(cor_unique, cor_gene1, cor_gene2) |> 
  select(-Var1, -Var2) |> 
  filter(gene1 != gene2) |> 
  arrange(desc(Freq)) 

# write out the results
write.csv(cor_df, "/Users/Ewann/splicing_comparison/results/nygc/cryptic_tdp_corr.csv", row.names = FALSE)
write.csv(cor_gene, "/Users/Ewann/splicing_comparison/results/nygc/cryptic_tdp_corr_by_gene.csv", row.names = FALSE)
