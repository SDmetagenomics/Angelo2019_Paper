library(DESeq2)
library(ggplot2)

######## END INITALIZATION #########




####################
#################### Differential Abundace Across Depth
####################




#### Load Data -- This will be the only part different between scripts analyzing different depths vs time/treatment ##
# setwd("~/Dropbox/Banfield_Lab_Files/Projects/Angelo2014/Studies/17_03_14_DEseq_of_m2_Mapping/")
# 
# countData <- read.delim("~/Dropbox/Banfield_Lab_Files/Projects/Angelo2014/Studies/17_03_13_Coverage_Results_from_m2_Mapping/Counts_Absolute_pruned.txt", header = FALSE,
#                         sep = "\t",row.names = 1) 
# 
# colData  <- read.delim("~/Dropbox/Banfield_Lab_Files/Projects/Angelo2014/Studies/17_03_14_DEseq_of_m2_Mapping/Sample_Metadata_ALL.txt", header = TRUE,
#                        sep = "\t") 

countData <- read.delim("Counts_Absolute_pruned.txt", header = FALSE, sep = "\t",row.names = 1) 

colData  <- read.delim("Sample_Metadata_ALL.txt", header = TRUE, sep = "\t") 



### Prepare Data

# Load Count and Metadata into DEseq Variables
countData <- as.matrix(countData)
colnames(countData) <- colData$Sample_Name # Use col1 as rownames
colData <- as.data.frame(colData)

# Make Metadata into Factors
colData$Time_Point <- as.factor(colData$Time_Point)
colData$Treat_Control <- as.factor(colData$Treat_Control)
colData$Replicate <- as.factor(colData$Replicate)
colData$Depth <- as.factor(colData$Depth)
colData$Plot <- as.factor(colData$Plot)




### Build DEseq Object - Depth
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ Replicate + Time_Point + Treat_Control + Depth)
print(dds)




### Run DEseq Stats with LRT
dds_lrt_out <- DESeq(dds, test = "LRT", 
                     reduced = ~ Replicate + Time_Point + Treat_Control)




### Output Size Factors and Plot Dispersion Estimates

# Output Size Factors
s_factors_lrt <- sizeFactors(dds_lrt_out)
#write.table(s_factors_lrt, "Depth_size_factors_lrt.txt", sep = "\t",quote = FALSE)

# Plot Dispersion Estimates
plotDispEsts(dds_lrt_out)




### Output Primary Results Table

# Show Results Variables
resultsNames(dds_lrt_out)

# Display Differential Abundance Summary
summary(results(dds_lrt_out, alpha = 0.05))

# Write Filtered Results (FDR < 0.05) To Disk
Depth_LRT_Results_p <- subset(results(dds_lrt_out, alpha = 0.05, format = "DataFrame"), padj < 0.05)
write.table(as.data.frame(Depth_LRT_Results_p), "Depth_LRT_Results_p.txt", sep = "\t", quote = FALSE)




### Run Linear Model on All Significant Results to Check Slope Across Depth

# Convert Results to DataFrame and add Colums for Slopes 
Results_Table <- as.data.frame(Depth_LRT_Results_p)
Results_Table <- data.frame(Results_Table, lm_slope = NA, slope_p = NA)

# Loop Through all Significant Results 
for (i in 1:length(rownames(Results_Table))){
  
  # Get Rowname of Row i
  tmp_OTU <- rownames(Results_Table)[i]
  
  # Get Counts vs Depth Data for SG, Log Scale Counts, and Create Numeric Dummy Depth Variable
  tmp_counts <- plotCounts(dds, gene=tmp_OTU, intgroup=c("Depth"),
                           returnData=TRUE)
  tmp_counts <- data.frame(tmp_counts,
                           log_count = log(tmp_counts$count), 
                           dummy_depth = as.numeric(tmp_counts$Depth))
  
  # Run Linear Model log_count = m*dummy_depth + b
  lm_d_log_dum <- lm(log_count ~ dummy_depth, data = tmp_counts)
  tmp_summary <- summary(lm_d_log_dum)
  
  # Add Model Coefficients to DataFrame
  Results_Table[i,7] <- tmp_summary$coefficients[2,1]
  Results_Table[i,8] <- tmp_summary$coefficients[2,4]
  
}

# FDR Correct Slope p-values and Subset for FDR <= 0.05
Results_Table <- data.frame(Results_Table, slope_fdr = p.adjust(Results_Table$slope_p, method = "fdr"))
Results_Table <- subset(Results_Table, slope_fdr <= 0.05)

# Add SG Centroid Names to Results table and Increase vs Decrease w/Depth Call
Results_Table <- data.frame(rpS3_Centroid = rownames(Results_Table),
                            Results_Table,
                            depth_call = ifelse(Results_Table$lm_slope > 0, "Increase", "Decrease"))





####################
#################### Differential Abundace Between Treatment v Control
####################




#### Load Data -- Specific Combined Variables for Compairing Treatment at Each Depth
#colData  <- read.delim("~/Dropbox/Banfield_Lab_Files/Projects/Angelo2014/Studies/17_03_14_DEseq_of_m2_Mapping/Sample_Metadata_ALL_Grouped.txt",
#                       header = TRUE, sep = "\t") 

colData  <- read.delim("Sample_Metadata_ALL_Grouped.txt", header = TRUE, sep = "\t") 



### Prepare Data

# Load New Metadata into DEseq Variables
colData <- as.data.frame(colData)

# Make Metadata into Factors
colData$Time_Point <- as.factor(colData$Time_Point)




### Build DEseq Object - Treatment 
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ Replicate + Time_Point + Factor)
dds



### Run DEseq Stats with Standard Wald Test with Local Fitting
dds_out <- DESeq(dds, fitType = "local")




### Output Size Factors and Plot Dispersion Estimates

# Output Size Factors
s_factors <- sizeFactors(dds_out)
#write.table(s_factors, "size_factors_local.txt", sep = "\t",quote = FALSE)

# Plot Dispersion Estimates
plotDispEsts(dds_out)




### Output Results Tables

# Show Results Variables
resultsNames(dds_out)

# Display Differential Abundance Summaries
summary(results(dds_out, contrast = c("Factor", "control20cm", "treatment20cm"), alpha = 0.05))
summary(results(dds_out, contrast = c("Factor", "control30cm", "treatment30cm"), alpha = 0.05))
summary(results(dds_out, contrast = c("Factor", "control40cm", "treatment40cm"), alpha = 0.05))

# Write Filtered Results (FDR < 0.05) To Disk
TvC_20cm <- results(dds_out, contrast = c("Factor", "control20cm", "treatment20cm"), alpha = 0.05) 
TvC_20cm_p <- subset(TvC_20cm, padj < 0.05)
#write.table(as.data.frame(TvC_20cm_p), "TvC_20cm_DESeq_ALL_local.txt", sep = "\t", quote = FALSE)
TvC_30cm <- results(dds_out, contrast = c("Factor", "control30cm", "treatment30cm"), alpha = 0.05) 
TvC_30cm_p <- subset(TvC_30cm, padj < 0.05)
#write.table(as.data.frame(TvC_30cm_p), "TvC_30cm_DESeq_ALL_local.txt", sep = "\t", quote = FALSE)
TvC_40cm <- results(dds_out, contrast = c("Factor", "control40cm", "treatment40cm"), alpha = 0.05)
TvC_40cm_p <- subset(TvC_40cm, padj < 0.05)
#write.table(as.data.frame(TvC_40cm_p), "TvC_40cm_DESeq_ALL_local.txt", sep = "\t", quote = FALSE)

