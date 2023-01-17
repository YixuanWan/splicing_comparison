library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
my_clean_reader <- function(file){
    df = data.table::as.data.table(janitor::clean_names(data.table::fread(file)))
    return(df)
}
fix_column_names = function(df){
    ori_colnames = df %>% colnames()
    ori_colnames[12] = "baseline_PSI"
    ori_colnames[13] = "contrast_PSI"
    colnames(df) = ori_colnames
    
    return(df)
    
}
#REPLACE WITH YOUR DOWNLOAD FOLDER PATH
file_path =file.path(here::here(),'majiq_sy5y_curves')
suffix = "_annotated_junctions.csv"
estimate_files = list.files(file_path,
                            pattern = suffix,
                            full.names = TRUE)

#Apply a function to each element of a list or atomic vector
sh_curve = purrr::map(estimate_files,my_clean_reader)

#simplify: list -> vector
samp_ids = base::tolower(purrr::simplify(purrr::map(estimate_files, basename)))
samp_ids = gsub(suffix,"",samp_ids)
##add on the name as an additional column
sh_curve = purrr::map2(sh_curve, samp_ids, ~cbind(.x, comparison = .y))
sh_curve = purrr::map(sh_curve,fix_column_names)
sh_curve = data.table::rbindlist(sh_curve)

sh_curve %>% 
    select(gene_name,paste_into_igv_junction,junc_cat,baseline_PSI,contrast_PSI,comparison,strand)  %>% 
    mutate(is_cryptic = contrast_PSI > 0.1 & baseline_PSI < 0.05)  %>%
    group_by(comparison,is_cryptic) %>% 
    summarize(n_jun = n_distinct(paste_into_igv_junction)) %>% 
    ungroup() %>% 
    filter(is_cryptic == TRUE) %>% 
    mutate(comparison = case_when(comparison == "nodox-dox00125" ~ "+",
                                  comparison == "nodox-dox00187" ~ "++",
                                  comparison == "nodox-dox0021" ~ "+++",
                                  comparison == "nodox-dox0025" ~ "++++",
                                  comparison == "nodox-dox0075" ~ "+++++")) %>% 
    ggplot(aes(x = comparison,
               y = n_jun)) + 
    geom_col() +
    # ggpubr::theme_pubr(base_size = 16) 
    theme_minimal() +
    ylab("N cryptic junctions") + 
    xlab(element_blank()) +
    ggtitle('N cryptic junctions in SH-SY5Y curve')
 
sh_curve %>% 
    select(gene_name,paste_into_igv_junction,junc_cat,,gene_id,baseline_PSI,contrast_PSI,comparison,strand) %>% 
    mutate(is_cryptic = contrast_PSI > 0.1 & baseline_PSI < 0.05) %>% 
    group_by(comparison,is_cryptic) %>% 
    summarize(n_jun = n_distinct(gene_id)) %>% 
    ungroup() %>% 
    filter(is_cryptic == TRUE) %>% 
    mutate(comparison = case_when(comparison == "nodox-dox00125" ~ "+",
                                  comparison == "nodox-dox00187" ~ "++",
                                  comparison == "nodox-dox0021" ~ "+++",
                                  comparison == "nodox-dox0025" ~ "++++",
                                  comparison == "nodox-dox0075" ~ "+++++")) %>% 
    ggplot(aes(x = comparison,
               y = n_jun)) + 
    geom_col() + 
    theme_minimal() +
    ylab("N cryptic genes") + 
    xlab(element_blank())


#REPLACE WITH YOUR DOWNLOAD FILE
file_path = file.path(here::here(),'majiq_dz_curves')
suffix = "_annotated_junctions.csv"
estimate_files = list.files(file_path,
                            pattern = suffix,
                            full.names = TRUE)

dz_curve = purrr::map(estimate_files,my_clean_reader)
dz_curve_samp_ids = base::tolower(purrr::simplify(purrr::map(estimate_files, basename)))
dz_curve_samp_ids = gsub(suffix,"",dz_curve_samp_ids)
##add on the name as an additional column
dz_curve = purrr::map2(dz_curve, dz_curve_samp_ids, ~cbind(.x, comparison = .y))
dz_curve = purrr::map(dz_curve,fix_column_names)
dz_curve = data.table::rbindlist(dz_curve)

dz_curve %>% 
    select(gene_name,paste_into_igv_junction,junc_cat,baseline_PSI,contrast_PSI,comparison,strand)  %>% 
    mutate(is_cryptic = contrast_PSI > 0.1 & baseline_PSI < 0.05)  %>% 
    group_by(comparison,is_cryptic) %>% 
    summarize(n_jun = n_distinct(paste_into_igv_junction)) %>% 
    ungroup() %>% 
    filter(is_cryptic == TRUE) %>% 
    mutate(comparison = case_when(comparison == "control-tdp43kd002" ~ "1",
                                  comparison == "control-tdp43kd003" ~ "2",
                                  comparison == "control-tdp43kd004" ~ "3",
                                  comparison == "control-tdp43kd005" ~ "4",
                                  comparison == "control-tdp43kd0075" ~ "5",
                                  comparison == "control-tdp43kd01" ~ "6",
                                  comparison == "control-tdp43kd1" ~ "7")) %>%
    ggplot(aes(x = comparison,
               y = n_jun)) + 
    geom_col() +
    # ggpubr::theme_pubr(base_size = 16) 
    theme_minimal() +
    ylab("N cryptic junctions") + 
    xlab(element_blank()) +
    ggtitle('N cryptic junctions in DZ curve')

dz_curve %>% 
    select(gene_name,paste_into_igv_junction,junc_cat,,gene_id,baseline_PSI,contrast_PSI,comparison,strand) %>% 
    mutate(is_cryptic = contrast_PSI > 0.1 & baseline_PSI < 0.05) %>% 
    group_by(comparison,is_cryptic) %>% 
    summarize(n_jun = n_distinct(gene_id)) %>% 
    ungroup() %>% 
    filter(is_cryptic == TRUE) %>% 
    mutate(comparison = case_when(comparison == "control-tdp43kd002" ~ "1",
                                  comparison == "control-tdp43kd003" ~ "2",
                                  comparison == "control-tdp43kd004" ~ "3",
                                  comparison == "control-tdp43kd005" ~ "4",
                                  comparison == "control-tdp43kd0075" ~ "5",
                                  comparison == "control-tdp43kd01" ~ "6",
                                  comparison == "control-tdp43kd1" ~ "7")) %>% 
    ggplot(aes(x = comparison,
               y = n_jun)) + 
    geom_col() + 
    theme_minimal() +
    ylab("N cryptic genes") + 
    xlab(element_blank())



dz_curve %>% 
    select(gene_name,paste_into_igv_junction,junc_cat,,gene_id,baseline_PSI,contrast_PSI,comparison,strand) %>% 
    filter(contrast_PSI > 0.1 & baseline_PSI < 0.05) 
