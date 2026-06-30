library(tidyr)
library(dplyr)
library(stringr)

# load file
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
if (str_sub(filename, -3) == "csv") {
	df <- read.csv(filename)
} else {
	df <- read.csv(paste(filename, "_raw.csv", sep = ""))
}

# get allele frequencies
df <- df %>%
  rowwise() %>%
  mutate(allele_freq = mean(c_across(X0:X99)))

df <- subset(df, allele_freq!=1.0 & allele_freq!=0.01)

foc_freq <- df$allele_freq[df$site_pos == 5000]
fix <- ifelse(length(foc_freq)==0, "none",
	ifelse(foc_freq > 0.35, "fix", "low_freq"))

df$afreq_bin <- ifelse(df$allele_freq >= 0.2 & df$allele_freq <= 0.8,
                       "common",
                       ifelse(df$allele_freq <= 0.05 | df$allele_freq >= 0.95,
                              "rare",
                              "int"))

d_ld_calc <- function(loc_row1, loc_row2){
  # d = pab - pa * pb
  pab <- sum(loc_row1 == loc_row2 & loc_row1 == 1)/length(loc_row1)
  D <- pab - (mean(loc_row1) * mean(loc_row2))
  
  if (D > 0){
    denom <- min(c((mean(loc_row1)*(1-mean(loc_row2))), (mean(1-loc_row1)*(mean(loc_row2)))))
  } else {
    denom <- min(c((mean(loc_row1)*(mean(loc_row2))), (mean(1-loc_row1)*(mean(1-loc_row2)))))
  }
  D_prime <- D/denom
  return(c(D, D_prime))
}

r2_ld_calc <- function(loc_row1, loc_row2, D){
  r2 <- D^2/(mean(loc_row1)*(1-mean(loc_row1))*mean(loc_row2)*(1-mean(loc_row2)))
  return(r2)
}

# function to compute r2
compute_ld_matrix <- function(geno_df) {
  cor_long <- expand.grid(snp1 = as.numeric(rownames(geno_df)), snp2 = as.numeric(rownames(geno_df)))
  cor_long <- subset(cor_long, snp1 < snp2)
  cor_long$dist <- cor_long$snp2 - cor_long$snp1
  cor_long$D <- NA
  cor_long$r2 <- NA
  cor_long$D_prime <- NA
  
  for (i in 1:nrow(cor_long)){
    loc_row1 <- as.numeric(geno_df[as.character(cor_long$snp1[i]),])
    loc_row2 <- as.numeric(geno_df[as.character(cor_long$snp2[i]),])
    dstats <- d_ld_calc(loc_row1, loc_row2)
    cor_long$D[i] <- dstats[1]
    cor_long$D_prime[i] <- dstats[2]
    cor_long$r2[i] <- r2_ld_calc(loc_row1, loc_row2, dstats[1])
  }
  
  return(cor_long)
}

# Create an empty dataframe and iterate over every category to find r2 between all pairs within said category
mut_types <- unique(df$site_type)
afreqs <- unique(df$afreq_bin)

linkage_df <- data.frame(snp1 = NA, snp2 = NA, dist = NA, D = NA, r2 = NA, D_prime = NA, site_type = NA, afreq_bin = NA, snp1_freq = NA, snp2_freq = NA)
for (mut in mut_types){
  for (a in afreqs){
    df_sub <- subset(df, site_type == mut & afreq_bin == a)
    if (nrow(df_sub) < 2){
      next
    } else {
      geno_df <- df_sub[,which(colnames(df_sub) == "X0"):which(colnames(df_sub) == "X99")]
      rownames(geno_df) <- df_sub$site_pos
      
      cor_long <- compute_ld_matrix(geno_df)
      cor_long$site_type <- mut
      cor_long$afreq_bin <- a
      cor_long2 <- merge(cor_long, df_sub[,c("site_pos", "allele_freq")], by.x = "snp1", by.y = "site_pos")
      cor_long2 <- merge(cor_long2, df_sub[,c("site_pos", "allele_freq")], by.x = "snp2", by.y = "site_pos")
      names(cor_long2)[names(cor_long2) == 'allele_freq.x'] <- 'snp1_freq'      
      names(cor_long2)[names(cor_long2) == 'allele_freq.y'] <- 'snp2_freq'
      
      linkage_df <- rbind(linkage_df, cor_long2)
    }
  }
}

# get rid of any NA values
linkage_df <- na.omit(linkage_df)
linkage_df$fix <- fix

# write output to csv
write.csv(linkage_df, paste(filename, "_pairwise", ".csv", sep = ""))
