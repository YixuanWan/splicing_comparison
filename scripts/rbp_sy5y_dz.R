library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(ggrepel)
library(tximport)
library(ggpubr)

cluster_mount_point = "/Users/Ewann/cluster"
tx2gene = file.path(cluster_mount_point,"vyplab_reference_genomes/annotation/human/GRCh38/gencode.v40.tx2gene.csv")
tx2gene <- data.table::fread(tx2gene,header=FALSE)  
colnames(tx2gene) = c("TXNAME", "GENEID")

human_rbp = data.table::fread("/Users/Ewann/splicing_comparison/human_rbps.csv",header=TRUE)


transcript_counts <- function(salmonpath,
                              metadatapath,
                              baseline = "control",
                              contrast = "contrast",
                              controls_name = "control",
                              contrast_name = "contrast",
                              countsFromAbundance = "lengthScaledTPM",
                              celltype) 
  {
  salmon_quant_directory = file.path(cluster_mount_point, salmonpath)
  metadata_filepath = file.path(cluster_mount_point, metadatapath)
  
  metadata_orig = read.csv(metadata_filepath,header=TRUE)
  column_name = "group"
  
  outputdir = "/Users/Ewann/splicing_comparison/rbp"
  exp = paste0(contrast_name,"_",controls_name)
  output_path = paste0(outputdir,"/",exp)
  
  metadata = metadata_orig %>%  
    filter(is.na(exclude_sample_downstream_analysis)) %>% 
    dplyr::select(sample_name, unit, !!(column_name)) %>%
    mutate(comparison_condition = case_when(!!as.symbol(column_name) == baseline ~ 'baseline',
                                            !!as.symbol(column_name) == contrast ~ 'contrast',
                                            TRUE ~ NA_character_)) %>% 
    filter(!is.na(comparison_condition))
  
  if(celltype == 'sy5y'){
    files = unique(file.path(salmon_quant_directory,paste0(metadata$unit,"_", metadata$sample_name),"quant.sf"))
  } else{
    files = unique(file.path(salmon_quant_directory,metadata$sample_name,"quant.sf")) 
  }
  
  names(files) = unique(metadata$sample_name)
  
  if(all(file.exists(files)) == FALSE) {
    stop("It seems that I cannot find those files...Please check if your directory is correct.")
  }
  
  txi.tx <- tximport(files, 
                     type="salmon", 
                     tx2gene=tx2gene,
                     ignoreTxVersion = TRUE,
                     ignoreAfterBar = TRUE,
                     txOut = TRUE)
  
  txi.sum <- summarizeToGene(txi.tx, tx2gene, countsFromAbundance = countsFromAbundance)
  genecount = txi.sum$abundance %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column('gene') %>% 
    mutate(ensmbleID = gsub("\\..*", "", gene)) %>% 
    left_join(human_rbp) %>% 
    filter(!is.na(gene_name))
  
  
  fwrite(genecount,paste0(output_path,".gene_counts.csv"))
}


#=========1. Generating transcript counts for DZ ========================== 

transcript_counts("first_weeks/TDP_curves_UPF1_GLIA/DZ_curves/salmon_quant",
                  "first_weeks/TDP_curves_UPF1_GLIA/DZ_curves/sample_sheet_dz_curves.csv",
                  baseline = "DZ_curves_0",
                  contrast = "DZ_curves_1",
                  controls_name = "control",
                  contrast_name = "DZhighest",
                  countsFromAbundance = "lengthScaledTPM",
                  "dz") 

#================2. Generating transcript counts for SY5Y =======================
transcript_counts("first_weeks/TDP_CHX_CLONES_GLIA/salmon",
                  "first_weeks/TDP_CHX_CLONES_GLIA/sample_sheet.csv",
                  baseline = "no_dox",
                  contrast = "dox_075",
                  controls_name = "control",
                  contrast_name = "SY5Yhighest",
                  countsFromAbundance = "lengthScaledTPM",
                  "sy5y") 

#================3. Identify and compare RBP expression in SY5Y and DZ =======================

sy5y = data.table::fread('/Users/Ewann/splicing_comparison/rbp/SY5Yhighest_control.gene_counts.csv')
dz = data.table::fread('/Users/Ewann/splicing_comparison/rbp/DZhighest_control.gene_counts.csv')


#======plotting top 10 expressed genes in sy5y ===============
sy5y_long = sy5y %>% 
  select(-gene) %>% 
  janitor::clean_names() %>% 
  melt(id.vars = c("gene_name","ensmble_id"))

avg_rbp_expression = sy5y_long %>% 
  group_by(gene_name) %>% 
  summarize(avg_expression = mean(value))
  
top10 = avg_rbp_expression %>% slice_max(avg_expression, n = 10) %>% pull(gene_name)


sy5y_long %>% 
  filter(gene_name %in% top10) %>%
  mutate(cond = ifelse(grepl("_nt_",variable),'control','TDPKD')) %>% 
  group_by(gene_name) %>% 
  mutate(avg_expresion = mean(value)) %>% 
  ungroup()  %>% 
  mutate(gene_name = forcats::fct_reorder(gene_name, avg_expresion)) %>% 
  ggplot(aes(y = value, x = gene_name)) +
  geom_boxplot() + 
  geom_jitter(height = 0,size = 3,aes(color = cond)) + 
  coord_flip() +
  ylab("Average expression of RBP in SH-SY5Y") +
  xlab("Gene name")+
  theme_minimal()

#======plotting top 10 expressed genes in dz ===========

dz_long = dz %>% 
  select(-gene) %>% 
  janitor::clean_names() %>% 
  melt(id.vars = c("gene_name","ensmble_id"))

avg_rbp_expression_d = dz_long %>% 
  group_by(gene_name) %>% 
  summarize(avg_expression = mean(value))

top10 = avg_rbp_expression_d %>% slice_max(avg_expression, n = 10) %>% pull(gene_name)


dz_long %>% 
  filter(gene_name %in% top10) %>%
  mutate(cond = ifelse(grepl("_0_",variable),'control','TDPKD')) %>% 
  group_by(gene_name) %>% 
  mutate(avg_expresion = mean(value)) %>% 
  ungroup()  %>% 
  mutate(gene_name = forcats::fct_reorder(gene_name, avg_expresion)) %>% 
  ggplot(aes(y = value, x = gene_name)) +
  geom_boxplot() + 
  geom_jitter(height = 0,size = 3,aes(color = cond)) + 
  coord_flip() +
  ylab("Average RBP expression in DZ") +
  xlab("Gene name") +
  theme_minimal()
  
#=========dz vs sy5y: overall ==========
rbpComp = left_join(avg_rbp_expression, avg_rbp_expression_d, by = "gene_name", suffix = c('.sy5y', '.dz')) %>% 
  as.data.frame() 
rbpComp %>% 
  ggplot(aes(x = avg_expression.dz, y = avg_expression.sy5y, label = gene_name)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm", size = 0.3) +
  scale_y_continuous(trans = scales::pseudo_log_trans()) +
  scale_x_continuous(trans = scales::pseudo_log_trans()) +
  geom_text_repel(box.padding = 0.2, max.overlaps = Inf, size = 2) +
  labs(title = 'Average RBP expression: DZ vs SH-SY5Y') +
  ylab("Average RBP expression/SH-SY5Y") +
  xlab("Average RBP expression/DZ") +
  theme_minimal() +
  geom_abline() +
  ggpubr::stat_cor()
  
#==========dz vs sy5y by condition =========

dz_tdpkd = dz_long %>% 
  mutate(cond = ifelse(grepl("_0_",variable),'control','TDPKD')) %>% 
  group_by(gene_name,cond) %>% 
  summarize(avg_expression = mean(value)) %>% 
  ungroup()

sy5y_tdpkd = sy5y_long %>% 
  mutate(cond = ifelse(grepl("_nt_",variable),'control','TDPKD')) %>% 
  group_by(gene_name,cond) %>% 
  summarize(avg_expression = mean(value)) %>% 
  ungroup()


rbp_tdpkd = left_join(sy5y_tdpkd, dz_tdpkd, by = c("gene_name","cond"), suffix = c('.sy5y', '.dz')) %>% 
  as.data.frame() 
p <- rbp_tdpkd %>% 
  ggplot(aes(x = avg_expression.dz, y = avg_expression.sy5y, label = gene_name)) +
  scale_y_continuous(trans = scales::pseudo_log_trans()) +
  scale_x_continuous(trans = scales::pseudo_log_trans()) +
  # geom_text_repel(box.padding = 0.2, max.overlaps = Inf, size = 2) +
  theme_minimal() +
  ylab("Average RBP expression/SH-SY5Y") +
  xlab("Average RBP expression/DZ") +
  geom_abline() +
  ggpubr::stat_cor()

p0 = p + 
  geom_point(size = 1, aes(color = cond)) +
  geom_smooth(method = "lm") +
  labs(title = 'DZ vs SH-SY5Y') 

p_con = p +
  geom_point(data = rbp_tdpkd[rbp_tdpkd$cond == "control",], size = 1) +
  geom_smooth(data = rbp_tdpkd[rbp_tdpkd$cond == "control",], method = "lm") +
  labs(title = 'DZ vs SH-SY5Y - control') 
  
p_kd = p +
  geom_point(data = rbp_tdpkd[rbp_tdpkd$cond == "TDPKD",], size = 1) +
  geom_smooth(data = rbp_tdpkd[rbp_tdpkd$cond == "TDPKD",], method = "lm") +
  labs(title = 'DZ vs SH-SY5Y - tdpkd') 

gridExtra::grid.arrange(p0, p_con, p_kd, ncol = 3)



