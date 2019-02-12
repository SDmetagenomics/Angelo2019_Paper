library(reshape)
library(plyr)
library(ggplot2)

#setwd("~/Dropbox/Banfield_Lab_Files/Manuscripts/2018_Angelo_C-N_Overview/FOR_SUBMISSION/Nat_Micro/Resubmission_1/Reviewer_Response/New_Analysis/Enrichment_Permutation_Test/")


### LOAD DATA ###
Phy_and_Targeted <- read.table("Phy_and_Targeted_Dat_S10.txt", header = T, sep = "\t", stringsAsFactors = F)



####
#### BEGIN FUNCTIONS ####
####

### 1
### Build Table for Phylum Testing
### removeZero removes functions/phyla that show zero frequency in the increased and decreased sets
build.phy.table <- function(x, fac, removeZero=T){
  tmp_phy_totals <- as.data.frame(table(x))
  tmp_phy_totals <- cast(tmp_phy_totals, formula(paste0("Phylum"," ~ ",fac)))
  tmp_phy_totals <- data.frame(tmp_phy_totals, 
                               In_Dec = sum(tmp_phy_totals$Decrease),
                               In_Inc = sum(tmp_phy_totals$Increase),
                               In_Nei = sum(tmp_phy_totals$Neither),
                               Total = apply(tmp_phy_totals, 1, sum))
  if (removeZero == T){
    test_zero <- tmp_phy_totals[,c(2,3)] # Test and Remove rows where Depth and Treatment Both have Zero
    test_zero <- which(rowSums(test_zero) == 0) # Test and Remove rows where Depth and Treatment Both have Zero
    
    if (length(test_zero) == 0){
      tmp_phy_totals <- tmp_phy_totals
    } else{
      tmp_phy_totals <- tmp_phy_totals[-test_zero,] # Test and Remove rows where Depth and Treatment Both have Zero
    }   
  }
  
  tmp_phy_totals
}


### 2
### Build Table for Targeted Function Testing
### removeZero removes functions/phyla that show zero frequency in the increased and decreased sets
build.targ.table <- function(x, fac, removeZero = T) {
  tmp_targ_totals <- aggregate(formula(paste0("."," ~ ",fac)), x, sum)
  tmp_targ_totals <- tmp_targ_totals[,-1]
  tmp_targ_totals <- t(tmp_targ_totals)
  tmp_targ_totals <- as.data.frame(tmp_targ_totals)
  tmp_targ_totals <- data.frame(gene = rownames(tmp_targ_totals),
                                Decrease = tmp_targ_totals$V1,
                                Increase = tmp_targ_totals$V2,
                                Neither = tmp_targ_totals$V3,
                                In_Dec = tmp_phy_totals[1,5],
                                In_Inc = tmp_phy_totals[1,6],
                                In_Nei = tmp_phy_totals[1,7],
                                Total = apply(tmp_targ_totals, 1, sum))
  rownames(tmp_targ_totals) <- NULL
  
  if (removeZero == T){
    test_zero <- tmp_targ_totals[,c(2,3)] # Test and Remove rows where Depth and Treatment Both have Zero
    test_zero <- which(rowSums(test_zero) == 0) # Test and Remove rows where Depth and Treatment Both have Zero
    
    if (length(test_zero) == 0){
      tmp_targ_totals <- tmp_targ_totals
    } else{
      tmp_targ_totals <- tmp_targ_totals[-test_zero,] # Test and Remove rows where Depth and Treatment Both have Zero
    }   
  }
  
  tmp_targ_totals
}


### 3
### Run Fisher Test with Post Hoc Hypergeometric Enrichment Test
test.fish <- function(x){
  test_table <- x
  test_table_out <- data.frame(test_table, Fish_All_p = NA, Log2_OR = NA)
  
  for (i in 1:nrow(test_table)){
    
    fish2x3_dat <- c(test_table[i,2], test_table[i,3], test_table[i,4], test_table[i,5] - test_table[i,2],
                     test_table[i,6] - test_table[i,3], test_table[i,7] - test_table[i,4])
    fish2x2_dat <- c(test_table[i,2], test_table[i,3], test_table[i,5] - test_table[i,2],
                     test_table[i,6] - test_table[i,3])
    
    fish2x3_dat <- matrix(fish2x3_dat, nrow = 2, ncol = 3, byrow = T,
                          dimnames = list(c("Yes", "No"), c("Decrease", "Increase", "Neither")))
    fish2x2_dat <- matrix(fish2x2_dat, nrow = 2, ncol = 2, byrow = T,
                          dimnames = list(c("Yes", "No"), c("Decrease", "Increase")))
    
    tmp_test_all <- fisher.test(fish2x3_dat,
                                alternative="two.sided")
    tmp_test_pair <- fisher.test(fish2x2_dat,
                                 alternative="two.sided")
    
    test_table_out[i,9] <- tmp_test_all$p.value
    test_table_out[i,10] <- log2(tmp_test_pair$estimate)
    
  }
  
  test_table_out <- data.frame(test_table_out,
                               Fish_All_FDR = p.adjust(test_table_out$Fish_All_p,
                                                       method = "fdr"))
  test_table_out
  
}


### 4
### Identify Comparisons for Permutation Test after Fischer Test and Run perm.test
perm.post.hoc <- function(tab, dat, group, nrep = 1000, type = c("binary", "category")){
  
  tab_filt <- subset(tab, Fish_All_FDR <= 0.1)
  
  perm_out <- data.frame(x = tab_filt[,1], abs_diff = NA, perm_p = NA)
  colnames(perm_out) <- c(colnames(tab_filt)[1], "abs_diff", "perm_p")
  
    
  if(type == "binary"){
    
    for(i in 1:nrow(perm_out)){
      
      perm_test_tmp <- perm.test(x = perm_out[i,1],
                                 dat = dat,
                                 group = group,
                                 nrep = nrep,
                                 type  = type)
      
      perm_out[i,2] <- perm_test_tmp[[2]]
      perm_out[i,3] <- perm_test_tmp[[5]]
      print(paste0("Completed Testing For: ",perm_out[i,1]))
    }
    
    perm_out <- data.frame(perm_out, perm_fdr = p.adjust(perm_out$perm_p))  
    
    merge_fac <- colnames(perm_out)[1]
    main_out <- merge(tab, perm_out, by = merge_fac, all.x = T)
  }
  
  
  if(type == "category"){
    
    category_col_name <- colnames(perm_out)[1]
    
    for(i in 1:nrow(perm_out)){
      
      perm_test_tmp <- perm.test(x = category_col_name,
                                 dat = dat,
                                 group = group,
                                 nrep = nrep,
                                 type  = type,
                                 cat_test = perm_out[i,1])
      
      perm_out[i,2] <- perm_test_tmp[[2]]
      perm_out[i,3] <- perm_test_tmp[[5]]
      print(paste0("Completed Testing For: ",perm_out[i,1]))
    }
    
    perm_out <- data.frame(perm_out, perm_fdr = p.adjust(perm_out$perm_p))  
    
    merge_fac <- colnames(perm_out)[1]
    main_out <- merge(tab, perm_out, by = merge_fac, all.x = T)
  }
  
  main_out

}


### 5
### Run Permutation test on selected function or Category
perm.test <- function(x, group, dat, nrep = 1000, type = c("binary", "category"), cat_test){
  
  
  ## Perform Testing for Binary Variable in Column
  if(type == "binary"){
  
    # Create List to hold output
    out_list <- list()
  
    # Parse Grouping and Functional Data Desired
    func_parse <- which(colnames(dat) == x)
    group_parse <- which(colnames(dat) == group)
    func_dat <- data.frame(group = dat[,group_parse],
                           x = dat[,func_parse])
    out_list[[1]] <- func_dat
  
    # Determine Native Counts and Groupings 
    sample_dist <- table(func_dat)
    group1 <- rownames(sample_dist)[1]
    group1_tot <- sample_dist[1,1] + sample_dist[1,2]
    group2 <- rownames(sample_dist)[2]
    group2_tot <- sample_dist[2,1] + sample_dist[2,2]
    
    # Calculate Original Fractional Difference for factor being tested
    group_hit_diff <- abs((sample_dist[1,2]/group1_tot) - (sample_dist[2,2]/group2_tot))
    
    grouping_vec <- c(rep(group1, group1_tot), rep(group2, group2_tot))
    
    # Perform Random Sampling and Group Assignment
    perm_res <- data.frame(Iteration = 1:nrep, Diff = 1:nrep)
    
    for(i in 1:nrep){
      perm_dat <- sample(func_dat$x,
                         size = length(grouping_vec),
                         replace = F)
      
      # deal with permutation that selects all zeros for function count
      if(sum(perm_dat) == 0){
        perm_res[i,2] <- 0
        next
      }
      
      group_dat <- sample(grouping_vec,
                          size = length(grouping_vec),
                          replace = F)
      test_dat <- data.frame(group = group_dat,
                             x = perm_dat)
      test_tab <- table(test_dat)
      
      # Calculate Test Statistic
      test_hit_diff <- ((test_tab[1,2]/group1_tot) - (test_tab[2,2]/group2_tot))
      
      # add permuted difference data to perm_res table
      perm_res[i,2] <- abs(test_hit_diff)
    }
    
    # Calculate p-value
    above_obs <- sum(perm_res$Diff >= group_hit_diff)
    p_val <- above_obs/nrep
    
    # Store All Data in Output List
    out_list[[2]] <- group_hit_diff
    out_list[[3]] <- perm_res
    out_list[[4]] <- above_obs
    out_list[[5]] <- p_val
  
  }
  
  ## Perform Testing for Categorical Variable in Column
  if(type == "category"){
    
    # Create List to hold output
    out_list <- list()
    
    # Parse Grouping and Functional Data Desired
    func_parse <- which(colnames(dat) == x)
    group_parse <- which(colnames(dat) == group)
    func_dat <- data.frame(group = dat[,group_parse],
                           x = dat[,func_parse])
    out_list[[1]] <- func_dat
    
    # Determine Native Counts and Groupings 
    sample_dist <- table(func_dat)
    group1 <- rownames(sample_dist)[1]
    group1_tot <- unname(rowSums(sample_dist)[1])
    group2 <- rownames(sample_dist)[2]
    group2_tot <- unname(rowSums(sample_dist)[2])
    
    group1_x_count <- sample_dist[1,which(colnames(sample_dist) == cat_test)]
    group2_x_count <- sample_dist[2,which(colnames(sample_dist) == cat_test)]
    
    # Calculate Original Fractional Difference for factor being tested
    group_hit_diff <- abs((group1_x_count/group1_tot) - (group2_x_count/group2_tot))
    
    grouping_vec <- c(rep(group1, group1_tot), rep(group2, group2_tot))
  
  
    # Perform Random Sampling and Group Assignment
    perm_res <- data.frame(Iteration = 1:nrep, Diff = 1:nrep)
  
    for(i in 1:nrep){
      perm_dat <- sample(func_dat$x,
                         size = length(grouping_vec),
                         replace = F)
      perm_dat <- ifelse(perm_dat == cat_test,1,0)
    
    # deal with permutation that selects all zeros for function count
      if(sum(perm_dat) == 0){
        perm_res[i,2] <- 0
        next
      }
    
      group_dat <- sample(grouping_vec,
                          size = length(grouping_vec),
                          replace = F)
      test_dat <- data.frame(group = group_dat,
                             x = perm_dat)
      test_tab <- table(test_dat)
      
      # Calculate Test Statistic
      test_hit_diff <- abs((test_tab[1,2]/group1_tot) - (test_tab[2,2]/group2_tot))
    
    # add permuted difference data to perm_res table
      perm_res[i,2] <- test_hit_diff
    }
    
    # Calculate p-value
    above_obs <- sum(perm_res$Diff >= group_hit_diff)
    p_val <- above_obs/nrep
    
    # Store All Data in Output List
    out_list[[2]] <- group_hit_diff
    out_list[[3]] <- perm_res
    out_list[[4]] <- above_obs
    out_list[[5]] <- p_val
  }
  
  out_list
  
}

####
#### END FUNCTIONS #### 
####





### Analyze Depth Gradient   

# Produce Count Tables Based on Depth
tmp_phy_totals <- build.phy.table(Phy_and_Targeted[,c(2,3)], "Depth", removeZero = T)
tmp_targ_totals <- build.targ.table(Phy_and_Targeted[,-c(1,2,4,5)], "Depth", removeZero = T) 

# Run Primary Fischer Testing
Depth_Targeted_Out <- test.fish(tmp_targ_totals)
Depth_Phy_Out <- test.fish(tmp_phy_totals)

# Run Post-Hoc Permutation Testing
set.seed(1234)
Depth_Targeted_Out <- perm.post.hoc(tab = Depth_Targeted_Out,
                                    dat = Phy_and_Targeted,
                                    group = "Depth",
                                    nrep = 10000,
                                    type = "binary")

set.seed(1234)
Depth_Phy_Out <- perm.post.hoc(tab = Depth_Phy_Out,
                               dat = Phy_and_Targeted,
                               group = "Depth",
                               nrep = 10000,
                               type = "category")



### Analyze Treatment Gradient - TvC 20cm

# Produce Count Tables Based on Depth
tmp_phy_totals <- build.phy.table(Phy_and_Targeted[,c(2,4)], "TvC20", removeZero = T)
tmp_targ_totals <- build.targ.table(Phy_and_Targeted[,-c(1,2,3,5)], "TvC20", removeZero = T) 

# Run Testing and Produce Final Output 
TvC20_Targeted_Out <- test.fish(tmp_targ_totals)
TvC20_Phy_Out <- test.fish(tmp_phy_totals)

# Run Post-Hoc Permutation Testing
set.seed(1234)
TvC20_Targeted_Out <- perm.post.hoc(tab = TvC20_Targeted_Out,
                                    dat = Phy_and_Targeted,
                                    group = "TvC20",
                                    nrep = 10000,
                                    type = "binary")

set.seed(1234)
TvC20_Phy_Out <- perm.post.hoc(tab = TvC20_Phy_Out,
                               dat = Phy_and_Targeted,
                               group = "TvC20",
                               nrep = 10000,
                               type = "category")



### Analyze Treatment Gradient - TvC 40cm

# Produce Count Tables Based on Depth
tmp_phy_totals <- build.phy.table(Phy_and_Targeted[,c(2,5)], "TvC40")
tmp_targ_totals <- build.targ.table(Phy_and_Targeted[,-c(1,2,3,4)], "TvC40") 

# Run Testing and Produce Final Output 
TvC40_Targeted_Out <- test.fish(tmp_targ_totals)
TvC40_Phy_Out <- test.fish(tmp_phy_totals)

# Run Post-Hoc Permutation Testing
set.seed(1234)
TvC40_Targeted_Out <- perm.post.hoc(tab = TvC40_Targeted_Out,
                                    dat = Phy_and_Targeted,
                                    group = "TvC40",
                                    nrep = 10000,
                                    type = "binary")

set.seed(1234)
TvC40_Phy_Out <- perm.post.hoc(tab = TvC40_Phy_Out,
                               dat = Phy_and_Targeted,
                               group = "TvC40",
                               nrep = 10000,
                               type = "category")



### Output Phylum Differential Abundance Tables
write.table(Depth_Phy_Out, "Depth_Phy_Out.txt", append = F, quote = F,row.names = F, sep = "\t")
write.table(Depth_Targeted_Out, "Depth_Targ_Out.txt", append = F, quote = F,row.names = F, sep = "\t")

write.table(TvC20_Phy_Out, "TvC20_Phy_Out.txt", append = F, quote = F,row.names = F, sep = "\t")
write.table(TvC20_Targeted_Out, "TvC20_Targ_Out.txt", append = F, quote = F,row.names = F, sep = "\t")

write.table(TvC40_Phy_Out, "TvC40_Phy_Out.txt", append = F, quote = F,row.names = F, sep = "\t")
write.table(TvC40_Targeted_Out, "TvC40_Targ_Out.txt", append = F, quote = F,row.names = F, sep = "\t")




### Produce Up Down Plots for Figure 4

phy_colors <- c(Acidobacteria = "#8DC63F", Actinobacteria = "#1B75BC",
                Alphaproteobacteria = "#939598", ANGP1 = "#BE1E2D",
                Bacteroidetes = "#00A79D", Bathyarchaeota = "#ED1C24",
                Betaproteobacteria = "#C2B59B", Chloroflexi = "#F16447",
                Deltaproteobacteria = "#8B5E3C", Euryarchaeota = "#F9ED32",
                Gammaproteobacteria = "#C49A6C", Gemmatimonadetes = "#EC008C",
                "NA" = "#000000", Nitrospirae = "#32C4E8",
                Planctomycetes = "#90AAE1", "RIF-WS3X" = "#006838",
                Rokubacteria = "#92278F", Thaumarchaeota = "#FBB040",
                Verrucomicrobia = "#39B54A")



### For Functions 
significant_functions <- c("URE","mauAB", "nit", "xoxF", "Frm2", "Fal1", "Fal2","coxL.I",
                           "amo_pmo", "nxrAB", "nirK")


# Subset Data for Significant Functions with Depth
Fig4B_Depth <- subset(Phy_and_Targeted, Depth == "Decrease" | Depth == "Increase")
Fig4B_Depth <- Fig4B_Depth[,-c(1,4,5)]
Fig4B_Depth <- aggregate(.~ Phylum + Depth, data = Fig4B_Depth, FUN = sum)
Fig4B_Depth <- melt(Fig4B_Depth)

Fig4B_Depth <- Fig4B_Depth[Fig4B_Depth$variable %in% significant_functions,]

depth_inc <- subset(Fig4B_Depth, Depth == "Increase")
depth_dec <- subset(Fig4B_Depth, Depth == "Decrease")

ggplot() +
  geom_bar(data = depth_dec, aes(x = variable, y = value, fill = Phylum),
           stat = "identity", color = "black", width = 0.8, size = 0.3) +
  geom_bar(data = depth_inc, aes(x = variable, y = -1*value, fill = Phylum),
           stat = "identity", color = "black",width = 0.8, size = 0.3) +
  geom_hline(yintercept = 0, size = 0.8) +
  scale_fill_manual(values = phy_colors) +
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position="none")


# Subset Data for Significant Functions with TvC20
Fig4B_TvC20 <- subset(Phy_and_Targeted, TvC20 == "Decrease" | TvC20 == "Increase")
Fig4B_TvC20 <- Fig4B_TvC20[,-c(1,3,5)]
Fig4B_TvC20 <- aggregate(.~ Phylum + TvC20, data = Fig4B_TvC20, FUN = sum)
Fig4B_TvC20 <- melt(Fig4B_TvC20)

Fig4B_TvC20 <- Fig4B_TvC20[Fig4B_TvC20$variable %in% significant_functions,]

TvC20_inc <- subset(Fig4B_TvC20, TvC20 == "Increase")
TvC20_dec <- subset(Fig4B_TvC20, TvC20 == "Decrease")

ggplot() +
  geom_bar(data = TvC20_dec, aes(x = variable, y = -1*value, fill = Phylum),
           stat = "identity", color = "black", width = 0.8, size = 0.3) +
  geom_bar(data = TvC20_inc, aes(x = variable, y = value, fill = Phylum),
           stat = "identity", color = "black",width = 0.8, size = 0.3) +
  geom_hline(yintercept = 0, size = 0.8) +
  scale_fill_manual(values = phy_colors) +
  xlab(NULL) +
  ylab(NULL) +
  #ylim(c(-125, 75)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position="none")


# Subset Data for Significant Functions with TvC40
Fig4B_TvC40 <- subset(Phy_and_Targeted, TvC40 == "Decrease" | TvC40 == "Increase")
Fig4B_TvC40 <- Fig4B_TvC40[,-c(1,3,4)]
Fig4B_TvC40 <- aggregate(.~ Phylum + TvC40, data = Fig4B_TvC40, FUN = sum)
Fig4B_TvC40 <- melt(Fig4B_TvC40)

Fig4B_TvC40 <- Fig4B_TvC40[Fig4B_TvC40$variable %in% significant_functions,]

TvC40_inc <- subset(Fig4B_TvC40, TvC40 == "Increase")
TvC40_dec <- subset(Fig4B_TvC40, TvC40 == "Decrease")

ggplot() +
  geom_bar(data = TvC40_dec, aes(x = variable, y = -1*value, fill = Phylum),
           stat = "identity", color = "black", width = 0.8, size = 0.3) +
  geom_bar(data = TvC40_inc, aes(x = variable, y = value, fill = Phylum),
           stat = "identity", color = "black",width = 0.8, size = 0.3) +
  geom_hline(yintercept = 0, size = 0.5) +
  scale_fill_manual(values = phy_colors) +
  xlab(NULL) +
  ylab(NULL) +
  ylim(c(-45, 20)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position="none")



### For Phyla
significant_phyla <- c("Gammaproteobacteria", "Euryarchaeota", "Bathyarchaeota", "Gemmatimonadetes",
                       "Rokubacteria", "Thaumarchaeota", "Chloroflexi", "Bacteroidetes",
                       "Verrucomicrobia")


# Subset Data for Significant Phyla with Depth
tmp_phy_totals <- build.phy.table(Phy_and_Targeted[,c(2,3)], "Depth", removeZero = F)
Fig4A_Dat <- tmp_phy_totals[tmp_phy_totals$Phylum %in% significant_phyla,]
Fig4A_Dat <- data.frame(Fig4A_Dat, Diff = (Fig4A_Dat$Decrease/Fig4A_Dat$In_Dec) - (Fig4A_Dat$Increase/Fig4A_Dat$In_Inc))
Fig4A_Dat <- Fig4A_Dat[,c(1,9)]

ordering <- data.frame(Fig4A_Dat[order(Fig4A_Dat$Diff, decreasing = T),], rank = 1:9)
ordering <- ordering[,-2]

ggplot(Fig4A_Dat, aes(x = reorder(Phylum, - Diff), y = Diff, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  scale_fill_manual(values = phy_colors) +
  geom_hline(yintercept = 0, size = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  ylim(c(-0.3,0.3)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position="none")


# Subset Data for Significant Phyla with TvC20
tmp_phy_totals <- build.phy.table(Phy_and_Targeted[,c(2,4)], "TvC20", removeZero = F)
Fig4A_Dat <- tmp_phy_totals[tmp_phy_totals$Phylum %in% significant_phyla,]
Fig4A_Dat <- data.frame(Fig4A_Dat, Diff = (Fig4A_Dat$Decrease/Fig4A_Dat$In_Dec) - (Fig4A_Dat$Increase/Fig4A_Dat$In_Inc))
Fig4A_Dat <- Fig4A_Dat[,c(1,9)]
Fig4A_Dat <- merge(Fig4A_Dat, ordering, by = "Phylum")


ggplot(Fig4A_Dat, aes(x = reorder(Phylum, rank), y = -1* Diff, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  scale_fill_manual(values = phy_colors) +
  geom_hline(yintercept = 0, size = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  ylim(c(-0.3,0.3)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position="none")


# Subset Data for Significant Phyla with TvC20
tmp_phy_totals <- build.phy.table(Phy_and_Targeted[,c(2,5)], "TvC40", removeZero = F)
Fig4A_Dat <- tmp_phy_totals[tmp_phy_totals$Phylum %in% significant_phyla,]
Fig4A_Dat <- data.frame(Fig4A_Dat, Diff = (Fig4A_Dat$Decrease/Fig4A_Dat$In_Dec) - (Fig4A_Dat$Increase/Fig4A_Dat$In_Inc))
Fig4A_Dat <- Fig4A_Dat[,c(1,9)]
Fig4A_Dat <- merge(Fig4A_Dat, ordering, by = "Phylum")

ggplot(Fig4A_Dat, aes(x = reorder(Phylum, rank), y = -1* Diff, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  scale_fill_manual(values = phy_colors) +
  geom_hline(yintercept = 0, size = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  ylim(c(-0.3,0.3)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position="none")


# Different Plot
# plot_dat <- subset(Depth_Targeted_Out, perm_fdr <= 0.05)
# plot_dat <- data.frame(plot_dat,
#                        Diff = ifelse(plot_dat$Log2_OR < 0, plot_dat$abs_diff * -1, plot_dat$abs_diff))
# 
# ggplot(plot_dat, aes(x = reorder(gene, - Diff), y = Diff)) +
#   geom_bar(stat = "identity", color = "black", width = 0.8, size = 0.3)


# ### FUNCTION TESTING ###
# x <- "Phylum"
# group <- "Depth"
# dat <- Phy_and_Targeted
# nrep <- 10000
# cat_test <- "Acidobacteria"
# tab <- Depth_Targeted_Out
# type <- "binary"
# 
# foo <- perm.test(x = "Phylum", group = "Depth", dat = Phy_and_Targeted, nrep = 1000, type = "category", cat_test = "Euryarchaeota")
# 
# foo <- perm.test(x = "nirK", group = "Depth", dat = Phy_and_Targeted, nrep = 1000, type = "binary")
# 
# foo
# ### FUNCTION TESTING ###

  