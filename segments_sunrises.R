# Produces the final segmentation results for HMF pipeline
# Written by Giulia 
# Sunrise plots and subclonality plots are created using https://github.com/hartwigmedical/hmftools/blob/master/purple/src/main/resources/r/somaticVariantPlots.R functions

####### INPUT:
####### 1. Sample_info_extended
####### 2. HOME
####### 3. publishDir

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(magrittr)
  library(ggplot2)
  library(RColorBrewer)
  library(ggfx)
})
theme_set(theme_bw())

conclusionLevels <- c("default", "adjust", "discard")
conclusionColors <- c("black", "#00e050", "#F04030")

#### Functions from hmftools repository
purity_ploidy_range_plot <- function(bestFit, range) {
  
  bestPurity = bestFit[1, "purity"]
  bestPloidy = bestFit[1, "ploidy"]
  bestScore = bestFit[1, "score"]
  
  range =  range %>%
    arrange(purity, ploidy) %>%
    group_by(purity) %>%
    mutate(
      absScore = pmin(4, score),
      score = pmin(1, abs(score - bestScore) / score),
      leftPloidy = lag(ploidy),
      rightPloidy = lead(ploidy),
      xmin = ploidy - (ploidy - leftPloidy) / 2,
      xmax = ploidy + (rightPloidy - ploidy) / 2,
      ymin = purity - 0.005,
      ymax = purity + 0.005,
      xmin = ifelse(is.na(xmin), ploidy, xmin),
      xmax = ifelse(is.na(xmax), ploidy, xmax))
  
  maxPloidy = min(range %>% arrange(purity, -ploidy) %>% group_by(purity)  %>% filter(row_number() == 1) %>% dplyr::select(purity, ploidy = xmax) %>% ungroup() %>% dplyr::select(ploidy))
  minPloidy = max(range %>% arrange(purity, ploidy) %>% group_by(purity)  %>% filter(row_number() == 1) %>% dplyr::select(purity, maxPloidy = xmin) %>% ungroup() %>% dplyr::select(maxPloidy))
  
  maxPloidy = max(maxPloidy, bestPloidy)
  minPloidy = min(minPloidy, bestPloidy)
  
  range = range %>%
    filter(xmin <= maxPloidy, xmax >= minPloidy) %>%
    mutate(xmax = pmin(xmax, maxPloidy), xmin = pmax(xmin, minPloidy))
  
  result = ggplot(range) +
    geom_rect(aes(fill=score, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
    scale_fill_gradientn(colours=c("blue","blue", "green", "yellow","orange", "red", "red2"), limits = c(0, 1), values=c(0, 0.0999, 0.1, 0.5, 0.8, 0.9, 1), breaks = c(0.1, 0.5, 1), labels = c("10%", "50%", "100%"), name = "Relative\nScore") +
    geom_segment(aes(y = 0.085, yend = 1.05, x=bestPloidy, xend = bestPloidy), linetype = "dashed", size = 0.1) +
    geom_label(data = data.frame(), aes(x = bestPloidy, y = 1.05, label = round(bestPloidy, 2)), size = 2.5) +
    geom_segment(aes(y = bestPurity, yend = bestPurity, x=minPloidy, xend = maxPloidy + 0.4), linetype = "dashed", size = 0.1) +
    geom_label(data = data.frame(), aes(y = bestPurity, x = maxPloidy + 0.4, label = paste0(bestPurity*100,"%" )), size = 2.5, hjust = 0.7) +
    theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position = "right", legend.title=element_text(size=6), legend.text=element_text(size=6)) +
    scale_y_continuous(labels = c("25%", "50%", "75%", "100%"), breaks = c(0.25, 0.5, 0.75, 1)) +
    xlab("Ploidy") + ylab("Purity") +  ggtitle("Purity/Ploidy Scores")
  
  return (result)
}

clonality_plot <- function(somaticBuckets, clonalityModel) {
  clonalityVariants = somaticBuckets %>% group_by(variantCopyNumberBucket) %>% summarise(count = sum(count))
  
  subclonalPercentage = clonalityModel %>% 
    group_by(bucket) %>% 
    mutate(totalWeight = sum(bucketWeight, na.rm = T)) %>%
    filter(isSubclonal) %>%
    summarise(
      isSubclonal = T,
      bucketWeight = sum(bucketWeight, na.rm = T), 
      subclonalLikelihood = ifelse(bucketWeight == 0 & is.nan(bucketWeight), 0, bucketWeight / max(totalWeight)))
  
  nonResidualModel = clonalityModel %>% filter(peak != 0)
  
  nonResidualSubclonalPercentage = nonResidualModel %>%
    group_by(bucket) %>%
    mutate(totalWeight = sum(bucketWeight)) %>%
    filter(isSubclonal) %>%
    summarise(
      isSubclonal = T,
      bucketWeight = sum(bucketWeight),
      subclonalLikelihood = ifelse(bucketWeight == 0, 0, bucketWeight / max(totalWeight)))
  
  combinedModel = nonResidualModel %>%
    group_by(bucket) %>% 
    summarise(bucketWeight = sum(bucketWeight))
  
  singleBlue = "#6baed6"
  singleRed = "#d94701"
  
  pTop = ggplot() +
    geom_bar(data=clonalityVariants, aes(x = variantCopyNumberBucket, weight = count), fill=singleBlue, col=singleBlue,  alpha = .4, size = 0.07, width = 0.05) +
    geom_line(data=combinedModel , aes(x = bucket, y = bucketWeight), position = "identity", alpha = 0.8) +
    geom_line(data=nonResidualModel, aes(x = bucket, y = bucketWeight, color = peak), position = "identity") +
    geom_area(data=nonResidualSubclonalPercentage %>% filter(isSubclonal), aes(x = bucket, y = bucketWeight), position = "identity",  alpha = 0.3, fill = singleRed, color = singleRed) +
    ggtitle("") + xlab("Variant Copy Number") + ylab("") +
    scale_y_continuous(expand=c(0.02, 0.02)) +
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position="none") +
    scale_x_continuous( expand=c(0.01, 0.01), limits = c(0, 3.5)) 
  
  pBottom = ggplot(data = subclonalPercentage) +
    geom_bar(width = 0.05, aes(x = bucket, y = subclonalLikelihood), stat = "identity", fill=singleRed, col=singleRed,  alpha = 0.3) + 
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank()) +
    xlab("") + ylab("") +
    #scale_y_continuous(labels = c("0%", "25%","50%","75%","100%"), breaks = c(0, 0.25, 0.5, 0.75, 1), expand=c(0.02, 0.02), limits = c(0, 1))# +
    scale_x_continuous( expand=c(0.01, 0.01), limits = c(0, 3.5)) 
  
  pFinal = cowplot::plot_grid(pTop, pBottom, ncol = 1, rel_heights = c(5, 1), align = "v")
  return(pFinal)
}

#### Candidates table
arg <- commandArgs(trailingOnly = TRUE)
sample_info<-suppressMessages(read_tsv(arg[1]))
publishDir <- arg[2]


# Selection of sample ids and patients of tumor samples 
sample_info <- sample_info %>%
  dplyr::filter(normal=="FALSE") %>%
  dplyr::select(sample, patient)

# Retrieve purity and ploidy values from purple results
sample_info$ploidy <- 0
sample_info$purity <- 0
sample_info$score <- 0
for (i in 1:nrow(sample_info)) {
  sample <- sample_info[[i, 1]]
  patient <- sample_info[[i, 2]]
  file <- paste0(publishDir, "/", patient, "/", sample, ".purple.purity.tsv")
  if (file.exists(file) == TRUE) {
    pur_pl <-read_tsv(file, show_col_types = FALSE)
    sample_info[i, 3] <- pur_pl[[1, 5]]
    sample_info[i, 4] <- pur_pl[[1, 1]]
    sample_info[i, 5] <- pur_pl[[1, 3]]
  }
}
# Check for completed samples
sample_info <- sample_info %>% mutate(sum = rowSums(.[3:5]))
samplesReady <- subset(sample_info, sum!=0)

# Join all segments retrieved from PURPLE results
#sample1 <- samplesReady[[1,1]]
#patient1 <- samplesReady[[1,2]]
all_segs <- data.frame()
  #read_tsv(paste0(HOME, patient, "/", patient, "_purple/", sample, ".purple.cnv.somatic.tsv"), show_col_types = FALSE)
#all_segs$sample <- sample1
#all_segs$patient <- patient1
for (i in 1:nrow(samplesReady)) {
  sample <- samplesReady[[i, 1]]
  patient <- samplesReady[[i, 2]]
  file <- paste0(publishDir, "/" , patient, "/", sample, ".purple.cnv.somatic.tsv")
  if (file.exists(file) == TRUE) {
    seg_part <- read_tsv(file, show_col_types = FALSE)
    seg_part$sample <- sample
    seg_part$patient <- patient
    all_segs <- rbind(all_segs, seg_part)
  }
}
# Calculare logR and LOH 
all_segs <- all_segs %>%
  inner_join(samplesReady) %>%
  mutate(length = end - start,
         logR = log(ifelse(copyNumber<0, 0, copyNumber)/ploidy),
         Loh = abs(baf - 0.5) * 2 * length / ifelse(is.na(baf), 0, length)) %>%
  dplyr::select(sample, patient, chromosome, start, end, length, 4:16, 19:20, 24:5)

write_tsv(all_segs, paste0(publishDir, "/segmentation_info_final.tsv"))

##### Sunrise and subclonal plots

all_candidates <- as.data.frame(matrix(ncol=5, nrow=0))
dir.create(paste0(publishDir, "/surnrises"))

for (i in 1:nrow(samplesReady)) {
  sample <- sample_info[i, 1]
  patient <- sample_info[i, 2]
  purpleDir <- paste0(publishDir, "/", patient)
  bestFitDF = read.table(file = paste0(purpleDir, "/", sample, ".purple.purity.tsv"), sep = "\t", header = T, comment.char = "!") %>% dplyr::select(purity, ploidy, score)
  rangeDF = read.table(file = paste0(purpleDir, "/", sample, ".purple.purity.range.tsv"), sep = "\t", header = T, comment.char = "!") %>%
    dplyr::select(purity, ploidy, score)
  rangePlot = purity_ploidy_range_plot(bestFitDF, rangeDF)
  ggsave(filename = paste0(sample, ".purity.range.png"), rangePlot, path = paste0(publishDir, "/sunrises/"), units = "in", height = 4, width = 4.8, scale = 1)
  
  somaticBuckets = read.table(paste0(purpleDir, "/", sample, ".purple.somatic.hist.tsv"), sep = "\t", header = T, numerals = "no.loss", skipNul = T)
  clonalityModel = read.table(paste0(purpleDir, "/", sample, ".purple.somatic.clonality.tsv"), sep = "\t", header = T, numerals = "no.loss", skipNul = T) %>%
    mutate(isSubclonal = isSubclonal == "true", isValid = isValid == "true", peak = as.character(peak), bucketWeight = as.numeric(as.character(bucketWeight))) %>% filter(isValid)
  clonalityModelPlot = clonality_plot(somaticBuckets, clonalityModel)
  ggsave(filename = paste0(sample, ".somatic.clonality.png"), clonalityModelPlot, path = paste0(publishDir, "/sunrises/"), units = "in", height = 6, width = 8, scale = 1)
}