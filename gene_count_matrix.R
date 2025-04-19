library(tidyr)
library(dplyr)
library(readr)

#add all csv files to list
files <- list.files('/home/mjw85/Documents/SRP/Group_Project/study_data',
                    pattern = '*.csv', full.name = TRUE)

files

#rename files in list by removing path info
names(files) <- stringr::str_split(files, pattern = '/', simplify = TRUE)[,8] %>%
  stringr::str_replace('.csv','')

files

#fill empty df with for loop
results<- data.frame()

for (file in files) {
  x<- read.csv(file, sep = '\t', header = FALSE, stringsAsFactors = FALSE)
  sample_name<- stringr::str_replace(file, '/home/mjw85/Documents/SRP/Group_Project/study_data', '') %>%
    stringr::str_replace('.csv','')
  x$sample<- sample_name
  results<- bind_rows(results, x)
}

results

#rename columns and pivot df to create gene count matrix
colnames(results)<- c('gene','expression','sample')

gene_matrix<- results %>%
  tidyr::pivot_wider(names_from = sample, values_from = expression)

gene_matrix

#create csv
write.csv(gene_matrix, '/home/mjw85/Documents/SRP/Group_Project/gene_count_matrix.csv')