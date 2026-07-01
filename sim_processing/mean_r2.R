# To run from command line
# Rscript --vanilla sim_processing/mean_r2_hitch.R $path/pairwise_ld_file.csv $S_BEN $S_DEL $RHO_EXP $MU_EXP $THETA $H_BEN $H_DEL $MODEL $NUM_SITES $output_file_name
# load libraries
library(dplyr)
library(tidyr)

# load in haplotype file
args <- commandArgs(trailingOnly = TRUE)
FILE_NAME <- args[1]
OUT_NAME <- args[11] #output file name should not include suffix (ie .csv) because we will make several output files that all start with the same root
df <- read.csv(FILE_NAME)

chr_length <- 10000

##### Step 1: get your classic mean r2 plot

# get fine scaled distance bin breaks.
breaks <- c(10, 15, 20, 25, 50, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500,
            550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1300, 1400, 1500, 1600, 1700, 1800,
            1900, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 7500, 10000)

mut_types <- c("nonsyn", "syn")
afreqs <- unique(df$afreq_bin) 

# construct empty output dataframe that yields all possible combinations of mutation type, allele frequency, and distance bin
mean_ld_df <- expand.grid(mut_types, afreqs, breaks)
colnames(mean_ld_df) <- c("site_type", "allele_freq", "dist_bin")
mean_ld_df$mean_D <- NA
mean_ld_df$mean_r2 <- NA
mean_ld_df$mean_D_prime <- NA
mean_ld_df$num_pairs <- NA

# calculate mean r2 for each combination of mutation type, allele frequency, and distance bin
for (i in 1:nrow(mean_ld_df)){
  # get distance minimum and maximum
  bin_max <- mean_ld_df$dist_bin[i]
  if (bin_max == min(breaks)) {
    bin_min <- 0
  } else {
    bin_min <- breaks[which(breaks == bin_max) - 1]
  }
  # take subset of dataframe
  df_sub <- subset(df, dist > bin_min & dist <= bin_max & site_type == mean_ld_df$site_type[i] & afreq_bin == mean_ld_df$allele_freq[i])
  # take average
  mean_ld_df$mean_D[i] <- mean(df_sub$D, na.omit = T)
  mean_ld_df$mean_r2[i] <- mean(df_sub$r2, na.omit = T)
  mean_ld_df$mean_D_prime[i] <- mean(df_sub$D_prime, na.omit = T)
  mean_ld_df$num_pairs[i] <- nrow(df_sub) # by incorporating number of pairs, if we want to use coarser bins later we can take a weighted average more easily
}

# add all parameter values to output file to make concatenating easier later
mean_ld_df$s_ben <- as.numeric(args[2]) 
mean_ld_df$s_del <- as.numeric(args[3])
mean_ld_df$rho_exp <- as.numeric(args[4])
mean_ld_df$mu_exp <- as.numeric(args[5])
mean_ld_df$theta <- as.numeric(args[6])
mean_ld_df$h_ben <- as.numeric(args[7])
mean_ld_df$h_del <- as.numeric(args[8])
mean_ld_df$model  <- as.numeric(args[9]) 
mean_ld_df$num_sites <- as.numeric(args[10])
write.csv(mean_ld_df, paste(OUT_NAME, "_mean_r2_global", ".csv", sep = ""))

##### Step 2: precise distance matching
paired <- df %>%
  group_by(dist, afreq_bin, snp1_freq, snp2_freq) %>%
  summarize(rN = mean(r2[site_type == "nonsyn"]),
            rS = mean(r2[site_type == "syn"]))

paired2 <- na.omit(paired)
paired2$delta_r2 <- paired2$rN - paired2$rS

paired2$s_ben <- as.numeric(args[2]) 
paired2$s_del <- as.numeric(args[3])
paired2$rho_exp <- as.numeric(args[4])
paired2$mu_exp <- as.numeric(args[5])
paired2$theta <- as.numeric(args[6])
paired2$h_ben <- as.numeric(args[7])
paired2$h_del <- as.numeric(args[8])
paired2$model  <- as.numeric(args[9]) 
paired2$num_sites <- as.numeric(args[10])
write.csv(paired2, paste(OUT_NAME, "_paired_delta_r2", ".csv", sep = ""))

##### Step 3: spatial LD in relationship to the focal allele
# add information about proximity to beneficial variant and whether both snps are on the same side of the beneficial variant or across it
df <- df %>% rowwise() %>%
  mutate(snp1_dist_focal = snp1 - 5000,
         snp2_dist_focal = snp2 - 5000,
		 min_abs_dist_focal = min(c(abs(snp1_dist_focal), abs(snp2_dist_focal))),
         max_abs_dist_focal = max(c(abs(snp1_dist_focal), abs(snp2_dist_focal))),
         across_ben = ifelse(snp1_dist_focal < 0 & snp2_dist_focal > 0, 
							 "across",
							 ifelse(snp1_dist_focal > 0 & snp2_dist_focal < 0, 
									"across", 
									"same")))

ben_dist_breaks <- c(15, 25, 50, 75)
ben_proximity <- c("near", "far")

spat_mean_ld_df <- expand.grid(mut_types, breaks, ben_dist_breaks, ben_proximity)
colnames(spat_mean_ld_df) <- c("site_type", "dist_bin", "ben_dist_breaks", "ben_proximity")
spat_mean_ld_df$mean_D <- NA
spat_mean_ld_df$mean_r2 <- NA
spat_mean_ld_df$mean_D_prime <- NA
spat_mean_ld_df$num_pairs <- NA

# calculate mean r2 for each combination of mutation type, distance bin, ben_dist_breaks, ben_proximity
# for this, we only care about afreq_bin == "common", across_ben == "same", and that the pair does not include the beneficial variant
for (i in 1:nrow(spat_mean_ld_df)){
  # get distance minimum and maximum
  bin_max <- spat_mean_ld_df$dist_bin[i]
  if (bin_max == min(breaks)) {
    bin_min <- 0
  } else {
    bin_min <- breaks[which(breaks == bin_max) - 1]
  }
  # take subset of dataframe
  if (spat_mean_ld_df$ben_proximity == "near"){
	  df_sub <- subset(df, dist > bin_min & dist <= bin_max & 
					   afreq_bin == "common" & across_ben == "same" & min_abs_dist_focal != 0 &
					   site_type == spat_mean_ld_df$site_type[i] & 
					   max_abs_dist_focal <= spat_mean_ld_df$ben_dist_breaks[i]) # for near we want both snps to be less than ben_dist_breaks away from beneficial variant
  } else {
	  df_sub <- subset(df, dist > bin_min & dist <= bin_max & 
					   afreq_bin == "common" & across_ben == "same" & min_abs_dist_focal != 0 &
					   site_type == spat_mean_ld_df$site_type[i] & 
					   min_abs_dist_focal > spat_mean_ld_df$ben_dist_breaks[i]) # for far we want both snps to be more than ben_dist_breaks away from beneficial variant
  }
	
  # take average
  spat_mean_ld_df$mean_D[i] <- mean(df_sub$D, na.omit = T)
  spat_mean_ld_df$mean_r2[i] <- mean(df_sub$r2, na.omit = T)
  spat_mean_ld_df$mean_D_prime[i] <- mean(df_sub$D_prime, na.omit = T)
  spat_mean_ld_df$num_pairs[i] <- nrow(df_sub) # by incorporating number of pairs, if we want to use coarser bins later we can take a weighted average more easily
}

spat_mean_ld_df$s_ben <- as.numeric(args[2]) 
spat_mean_ld_df$s_del <- as.numeric(args[3])
spat_mean_ld_df$rho_exp <- as.numeric(args[4])
spat_mean_ld_df$mu_exp <- as.numeric(args[5])
spat_mean_ld_df$theta <- as.numeric(args[6])
spat_mean_ld_df$h_ben <- as.numeric(args[7])
spat_mean_ld_df$h_del <- as.numeric(args[8])
spat_mean_ld_df$model  <- as.numeric(args[9]) 
spat_mean_ld_df$num_sites <- as.numeric(args[10])
write.csv(paired2, paste(OUT_NAME, "_spatial_ld", ".csv", sep = ""))
