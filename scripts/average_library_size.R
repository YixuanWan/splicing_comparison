library(dplyr)
library(tidyr)
my_clean_reader <- function(file){
  df = data.table::as.data.table(janitor::clean_names(data.table::fread(file)))
  return(df)
}

lib_mean <- function(file){
  result = mean(file$library_size)
  return(result)
}

experiment_table <- read.csv("TDP-43 KD experiments - Sheet1.csv")

#REPLACE WITH YOUR DOWNLOAD FILE
file_path = file.path(here::here(),'library_size')
prefix = "library_size_"
suffix = ".csv"
estimate_files = list.files(file_path,
                            pattern = prefix,
                            full.names = TRUE)

#cleaning the data
libsize = purrr::map(estimate_files,my_clean_reader)

#creating a list of mean values
average = purrr::map(libsize, lib_mean) %>% purrr::simplify()

#extracting the experiment names
experiment_names = base::tolower(purrr::simplify(purrr::map(estimate_files, basename)))
experiment_names = gsub(suffix,"",experiment_names)
experiment_names = gsub(prefix,"",experiment_names) 

#merging averages with experiment names
tibble(experiment_names,average)


  