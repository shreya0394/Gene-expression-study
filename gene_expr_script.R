# manipulating gene expression data

# Load the required libraries
library(tidyverse)
library(dplyr)
library(GEOquery)
library(ggplot2)
library(readxl)

# Reading input data file
data_file <- read_excel("GSE229339_Complete_Data_FPKM.xlsx")
head(data_file)

# Getting metadata
dat_gse <- getGEO(GEO = "GSE229339", GSEMatrix = TRUE)
head(dat_gse)

# Extracting phenodata
metadata <- pData(phenoData(dat_gse[[1]]))

# Modifying obtained metadata
modified_met <-metadata %>%
  select(1,11,12,13,19) %>%
  rename(cell_line = 'characteristics_ch1.1', treatment = 'characteristics_ch1.2', genotype = 'characteristics_ch1.3') %>%
  mutate(cell_line = gsub("cell line: ",'', cell_line), treatment= gsub("treatment: ",'', treatment),genotype = gsub("genotype: ",'',genotype))

# Reshaping data
long_dat<- data_file %>%
  rename(genes = 1) %>%
  gather(key = 'samples', value = 'FPKM', -genes)

colnames(modified_met)
colnames(long_dat)

# Joining data
result <- long_dat %>%
  left_join(.,modified_met, by = c('samples' = 'description'))

# Exploring data, check for any missing values
any(is.na(result))

# Finding outliers in data
mean_FPKM <- mean(result$FPKM, na.rm = TRUE)
sd_FPKM <- sd(result$FPKM, na.rm = TRUE)

# Calculate Z-scores for FPKM values to find outliers with threshold -3 to 3
Z_score <- (result$FPKM - mean_FPKM) / sd_FPKM
Z_score

outliers <- result[result$Z_score > 3 | result$Z_score < -3, ]
outliers

# Statistical analysis
# t-test for comparing treatment and FPKM values
group1 <- result[result$treatment == "UNTREATED", "FPKM"]
group2 <- result[result$treatment == "1.3 mer transfected", "FPKM"]

ttest_result <- t.test(group1, group2)
ttest_result

# Compare distribution of data based on treatment
wilcox.test(FPKM ~ treatment, data = result)


# Filtering data by required parameter
filtered_result <- result %>%
  filter(genotype == 'control' | genotype == 'UBR7 knockdown') %>%
  group_by(genes, treatment) %>%
  summarise(
    mean_FPKM = mean(FPKM),
    median_FPKM = median(FPKM)
  ) %>%
  arrange(-mean_FPKM)

# Data visualization
# bar plot
result %>%
  filter(genotype == 'UBR7 knockdown'| genotype == 'control') %>%
  ggplot(aes(x = samples, y = FPKM, fill = treatment)) +
  geom_col()

# density plot
result %>%
  filter(genes=='RPLP1') %>%
  ggplot(aes(x = FPKM, fill = genotype))+
  geom_density(alpha = 0.5)

# heatmap
# gene of interest
interest <- c('RNR2', 'GADPH', 'FTH1', 'RPLP1', 'RPL8', 'RPS11', 'ACTG1', 'ALDOA')

new_filtered_result <- result %>%
  filter(genes %in% interest)%>%
  head()

plot <- new_filtered_result %>%
  ggplot(aes(x = samples, y = genes, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'blue') +
  geom_text(aes(label = sprintf("%.1f", FPKM)), vjust = 1, size=3,alpha = 0.5) + 
  theme_minimal()

# Print the plot
print(plot)  



  
  
