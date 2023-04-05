library(arrow)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)


# load in dataset
nygc = read_parquet("/Users/Ewann/splicing_comparison/data/selective_cryptic_psi_in_nygc.parquet")

psi_tdp = nygc |>
  filter(tdp_path == "path") |> 
  select(psi, paste_into_igv_junction) |> 
  group_by(paste_into_igv_junction) |> 
  pivot_wider(names_from = paste_into_igv_junction, 
              names_repair = "check_unique", 
              values_from = psi,
              values_fn = list) |> 
  purrr::map(as.data.table)

psi_wide = do.call(cbind, psi_tdp) |> setnames(names(psi_tdp))

cor_matrix = psi_wide |> cor() 
diag(cor_matrix) = 0

# Convert the correlation matrix to a long format data frame
cor_df = as.data.frame(as.table(cor_matrix))|> 
  arrange(desc(Freq)) 

# Identify the most correlated cryptics
top_cor = cor_df |>  
  filter(Freq > 0.5 | Freq < -0.5)

# Create the heatmap using ggplot2
top_cor |>   
ggplot(aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "lightgrey", high = "purple") +
  labs(x = "", y = "", fill = "Correlation") +
  theme_minimal() +
  # theme(axis.text.x = element_blank(), axis.text.y = element_blank())
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# write out the results
write.csv(cor_df, "/Users/Ewann/splicing_comparison/results/nygc/cryptic_tdp_corr.csv", row.names = FALSE)
  