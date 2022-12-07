library(dplyr)
library(magrittr)
library(readr)
library(tibble)
library(tidyverse)

arg <- commandArgs(trailingOnly = TRUE)
#arg <- "/Users/micoli/mnt/storageBig8/resources/processed_data/HERCULES/WGSbams/sample_info_extended.csv"
#proj_dir <- "/Users/micoli/mnt/storageBig8/work/micoli/pipeline/"

samples <- read_tsv(arg[1], show_col_types = FALSE)
samples$normalSample<-replace_na(samples$normalSample, "unknown")
proj_dir <- arg[2]

# File for all tools
ids_table <- samples %>%
  filter(normal!=TRUE) %>%
  select(sample, normalSample, bamFile, normalBamFile, patient)

## File for Amber and Cobalt 
ids_table_filtered <- data.frame(matrix(ncol=ncol(ids_table), nrow = 0))
colnames(ids_table_filtered) <- colnames(ids_table)
for (i in 1:nrow(ids_table))
  {
    sample <- ids_table[[i, 1]]
    patient <- ids_table[[i, 5]]
    path <- paste0(proj_dir,"/", patient, "/", sample, ".cobalt.ratio.pcf")
    if (file.exists(path)==F)
    {
      ids_table_filtered <- rbind(ids_table[i, ], ids_table_filtered)
    }
  }
write_tsv(ids_table_filtered, "ids_table_filtered.tsv", na = "")

## File for Gripss and Purple
patients_redo <- unique(as.data.frame(ids_table_filtered[,5]))
ids_table_purple <- ids_table %>%
  filter(patient %in% patients_redo$patient)
write_tsv(ids_table_purple, "ids_table_gripss.tsv")

# Prepare file for GRIDSS: normal and all tumor bams related and put in line
# Select the non-processed samples
normal_samples <- samples[samples$normal==TRUE,]
normal_samples <- normal_samples %>%
  filter(patient %in% ids_table_filtered$patient)
gridss_input <- normal_samples[,c(1, 2, 13)]

gridss_input$tumorBam <- ""
for (i in gridss_input$sample) {
  j <- NULL
  j <- as.data.frame(cbind(j, ifelse(samples$normalSample == i & samples$normal == "FALSE", samples$bamFile, NA))) %>%
    na.omit(j)
  j<- toString(c(j$V1[1:nrow(j)]))
  j<- str_remove_all(j, ",")
  
  gridss_input$tumorBam <- ifelse(gridss_input$sample == i, j, gridss_input$tumorBam)
}

gridss_input <- gridss_input[order(gridss_input$sample),] 
gridss_input <- as.data.frame(gridss_input)
gridss_input$tumorBam <- str_replace(gridss_input$tumorBam, "^NA$", "")

#The samples which have a different normal (different platforms) but belong to the same patient are joined together to 
#increase the power of the joint calling
duplicato<- gridss_input[duplicated(gridss_input$patient),] #Other platform samples
true_normal <- gridss_input[!duplicated(gridss_input$patient),]  #Normal BDNA 
duplicato <- duplicato %>% unite("Bams", bamFile:tumorBam, sep=" ", remove = FALSE) %>% select(patient, Bams)

#If there are patients with multiple normal Bams they are joined to the true normal BDNA
if (nrow(duplicato > 0)) {
  true_normal$other_plat<-""
  for (line in 1:nrow(duplicato)) {
    pat <- duplicato[line, 1]
    bams <- duplicato[line, 2]
    row<-grep(pat, true_normal$patient)
    true_normal[row, "other_plat"]<-bams
  }
  
  true_normal$other_plat <- replace_na(true_normal$other_plat, "")
  true_normal <- true_normal %>% unite("all_bams", bamFile:other_plat, sep=" ", remove = FALSE) %>% select(sample, patient, all_bams)
  colnames(true_normal)<- c("normalSample", "patient", "Bams")
  write_tsv(true_normal, "gridss_input.tsv", na = "")
  
} else {
  gridss_input <- gridss_input %>% unite("all_bams", bamFile:tumorBam, sep=" ", remove = FALSE) %>% select(sample, patient, all_bams)
  colnames(gridss_input)<- c("normalSample", "patient", "Bams")
  write_tsv(gridss_input, "gridss_input.tsv", na = "")
}

