library(dplyr)
library(tidyr)

# load file
# Rscript --vanilla mean_r2_hitch.R $TARGET_FOLDER $SUBDIR_NAME $output_file $S_TOG $THETA $S_DEL $RHO_EXP $MUT_DEL $CHR_LENGTH
args <- commandArgs(trailingOnly = TRUE)
BIG_FOLDER <- args[1]
SUBDIR_NAME <- args[2]
FILE_NAME <- args[3]
s_focal <- as.numeric(args[4])
theta <- as.numeric(args[5])
model  <- as.numeric(args[6])
x0 <- as.numeric(args[7])
mean_ep_dist <- as.numeric(args[8])
h_ben <- as.numeric(args[9])
h_del <- as.numeric(args[10])
end_freq <- as.numeric(args[11])
chr_length <- as.numeric(args[12])
df <- read.csv(paste(BIG_FOLDER, "/", SUBDIR_NAME, "/", FILE_NAME, sep = ""))

##### Step 1 + 2: get your classic mean r2 plot. this is also going to decompose rN into pairs that include the beneficial allele and those that don't

# get breaks that are evenly spaced on a log scale
breaks <- c(10, 15, 20, 25, 50, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500,
            550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1300, 1400, 1500, 1600, 1700, 1800,
            1900, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 7500, 10000)

# list out all unique combos of mutation type and allele frequency
ben_loc <- chr_length/2 
df$site_type <- ifelse(df$snp1 == ben_loc | df$snp2 == ben_loc, "nonsyn_ben", df$site_type)
mut_types <- c("nonsyn", "nonsyn_ben", "nonsyn_all", "syn")
#mut_types <- unique(df$site_type)
afreqs <- unique(df$afreq_bin)
table(df$site_type)

# construct output dataframe that yields all possible combinations of mutation type, allele frequency, and distance bin
mean_ld_df <- expand.grid(mut_types, afreqs, breaks)
colnames(mean_ld_df) <- c("site_type", "allele_freq", "dist_bin")
mean_ld_df$mean_D <- NA
mean_ld_df$mean_r2 <- NA
mean_ld_df$mean_D_prime <- NA

# calculate mean r2 for each set of parameters
for (i in 1:nrow(mean_ld_df)){
  bin_max <- mean_ld_df$dist_bin[i]
  if (bin_max == min(breaks)) {
    bin_min <- 0
  } else {
    bin_min <- breaks[which(breaks == bin_max) - 1]
  }
  df_sub <- subset(df, dist > bin_min & dist <= bin_max & site_type == mean_ld_df$site_type[i] & afreq_bin == mean_ld_df$allele_freq[i])
  mean_ld_df$mean_D[i] <- mean(df_sub$D, na.omit = T)
  mean_ld_df$mean_r2[i] <- mean(df_sub$r2, na.omit = T)
  mean_ld_df$mean_D_prime[i] <- mean(df_sub$D_prime, na.omit = T)
  mean_ld_df$num_pairs[i] <- nrow(df_sub)
}

# extract values of parameters
mean_ld_df$s_focal <- s_focal
mean_ld_df$theta <- theta
mean_ld_df$s_del <- s_del
mean_ld_df$recomb_exp <- recomb_exp
mean_ld_df$mut_exp <- mut_del
mean_ld_df$h_ben <- h_ben
mean_ld_df$h_del <- h_del
mean_ld_df$end_freq <- end_freq

write.csv(mean_ld_df, paste(BIG_FOLDER, "/", SUBDIR_NAME, "_mean_r2_global", ".csv", sep = ""))

##### Step 3: precise distance matching
paired <- df %>%
  group_by(dist, afreq_bin, snp1_freq, snp2_freq) %>%
  summarize(rN = mean(r2[site_type == "nonsyn"]),
            rS = mean(r2[site_type == "syn"]))

paired2 <- na.omit(paired)
paired2$delta_r2 <- paired2$rN - paired2$rS

paired2$s_focal <- s_focal
paired2$theta <- theta
paired2$recomb_exp <- recomb_exp

write.csv(paired2, paste(BIG_FOLDER, "/", SUBDIR_NAME, "_paired_delta_r2", ".csv", sep = ""))

##### Step 4: spatial LD in relationship to the focal allele
df <- df %>% rowwise() %>%
  mutate(snp1_dist_focal = snp1 - ben_loc,
         snp2_dist_focal = snp2 - ben_loc,
	 snp_dist = abs(snp1-snp2),
         min_abs_dist_focal = min(c(abs(snp1_dist_focal), abs(snp2_dist_focal))),
         max_abs_dist_focal = max(c(abs(snp1_dist_focal), abs(snp2_dist_focal))),
         mean_dist_focal = mean(c(snp1_dist_focal, snp2_dist_focal)), 
         across_ben = ifelse(snp1_dist_focal < 0 & snp2_dist_focal > 0, "across",
                              ifelse(snp1_dist_focal > 0 & snp2_dist_focal < 0, "across", "same")))

write.csv(subset(df, df$afreq_bin  == "common"), paste(BIG_FOLDER, "/", SUBDIR_NAME, "_all_pairs", ".csv", sep = ""))

#spatial_mean_ld_df <- expand.grid(mut_types, afreqs, unique(df$across_ben), breaks[2:length(breaks)], c("near", "far"))
#colnames(spatial_mean_ld_df) <- c("site_type", "allele_freq", "across_ben", "dist", "dist_from_ben")
#spatial_mean_ld_df$mean_r2 <- NA
#spatial_mean_ld_df$num_pair <- NA
#
#for (i in 1:nrow(spatial_mean_ld_df)){
#  ben_bin_max <- spatial_mean_ld_df$dist_from_ben[i]
#  ben_bin_min <- spatial_breaks[which(spatial_breaks == ben_bin_max) - 1]

#  bin_max <- spatial_mean_ld_df$dist_bw[i]
#  if (bin_max == min(dist_breaks)) {
#    bin_min <- 0
#  } else {
#    bin_min <- dist_breaks[which(dist_breaks == bin_max) - 1]
#  }  
#
#  if (spatial_mean_ld_df$site_type[i] == "nonsyn_all") {
#    df_sub <- subset(df, mean_dist_focal > ben_bin_min & mean_dist_focal <= ben_bin_max & snp_dist > bin_min & snp_dist <= bin_max & site_type != "syn" & afreq_bin == spatial_mean_ld_df$allele_freq[i] & across_ben == spatial_mean_ld_df$across_ben[i])
#  } else {
#    df_sub <- subset(df, mean_dist_focal > ben_bin_min & mean_dist_focal <= ben_bin_max & snp_dist > bin_min & snp_dist <= bin_max & site_type == spatial_mean_ld_df$site_type[i] & afreq_bin == spatial_mean_ld_df$allele_freq[i] & across_ben == spatial_mean_ld_df$across_ben[i])
#  }
#
#  spatial_mean_ld_df$mean_r2[i] <- mean(df_sub$r2, na.omit = T)
#  spatial_mean_ld_df$num_pair[i] <- nrow(df_sub)
#}

# extract values of parameters
#spatial_mean_ld_df$s_focal <- s_focal
#spatial_mean_ld_df$theta <- theta
#spatial_mean_ld_df$recomb_exp <- recomb_exp

#write.csv(spatial_mean_ld_df, paste(BIG_FOLDER, "/", SUBDIR_NAME, "_spatial_ld", ".csv", sep = ""))
