library(dplyr)
library(tidyr)
library(VennDiagram)

splicingTable = read.csv("/Users/Ewann/splicing_comparison/data/majiq/splicing_full_delta_psi_tables.csv")
metatable = read.csv("/Users/Ewann/splicing_comparison/samplesheet/tdp_experiments_updated-tdp_experiments_updated.csv")


crypticTable = splicingTable |> 
  select(gene_name,paste_into_igv_junction,junc_cat,baseline_PSI,contrast_PSI,comparison,strand) |>  
  filter(comparison %in% metatable$comparison) |> 
  # filter(junc_cat != "annotated") |>
  mutate(is_cryptic = contrast_PSI > 0.1 & baseline_PSI < 0.05)  |> 
  filter(is_cryptic == TRUE) |> 
  group_by(paste_into_igv_junction, comparison) |> 
  reframe(gene_name, junc_cat, baseline_PSI = mean(baseline_PSI), contrast_PSI = mean(contrast_PSI), comparison, strand, is_cryptic) |> 
  ungroup() |> 
  distinct(paste_into_igv_junction, comparison, .keep_all = TRUE)

# find overlaps across cell types
celltype_overlap = crypticTable |> 
  filter(comparison %in% c("controlbrowncorticalneuron-tdp43kdbrowncorticalneuron", 
                           "controlhumphreycorticalneuron-tdp43kdhumphreycorticalneuron",
                           "controlseddighicorticalneuron-tdp43kdseddighicorticalneuron",
                           "nodox-dox0075", "control-tdp43kd1",
                           "controlklimmotorneuron-tdp43kdklimmotorneuron",
                           "controlbrownshsy5y-tdp43kdbrownshsy5y",
                           "controlbrownskndz-tdp43kdbrownskndz")) |> 
  select(comparison, paste_into_igv_junction) |> 
  unique() |> 
  pivot_wider(names_from = "comparison",
              values_from = "paste_into_igv_junction")

celltype_overlap = purrr::map(celltype_overlap, unlist)

# Create the Venn diagram
venn.diagram(celltype_overlap[c(3,4,8)], 
             filename = NULL,
             category.names = c("SH-SY5Y", "i3Neuron", "SK-N-DZ"),
             fill = c("#619CFF","salmon", "#E76BF3"),
             col = "lightgrey",
             cex = 0.8,
             cat.cex = 0.8,
             alpha = 0.5,
             main = "Cryptic Splicing Junctions") |> grid.draw()

