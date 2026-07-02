library(tidyverse)

##### Make iLDS scan figure (4A)
# Read in output from iLDS
ilds <- read.delim("C_diff_full_scan.txt", sep = ",")
# Calculate Signed iLDS
ilds$silds <- ilds$rNrS * ilds$r2L/begg_ilds$r2G

# Grey rectangles around two focal sweeps
sweep_centers <- c(2740747, 2976139)
dec_dist <- 10^4.5
rects_ilds <- data.frame(
  xmin = c(sweep_centers[1]-dec_dist, sweep_centers[2]-dec_dist),
  xmax = c(sweep_centers[1]+dec_dist, sweep_centers[2]+dec_dist),
  ymin = c(-0.25, -0.25),
  ymax = c(0.75, 0.75),
  fill_color = c("grey", "grey")
)

ildsplot <- ggplot(data = subset(ilds)) + 
  geom_rect(data = rects_ilds, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "lightgrey", alpha = 0.5) +
  geom_point(aes(x = site_pos, y = silds, color = significance)) + 
  scale_color_manual(values = c("False" = "lightblue", "True" = "orange")) +
  theme_bw() + 
  geom_hline(yintercept = 0) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 18)) +
  ylab("Signed iLDS (rN-rS) * rL/rG") +
  xlab("Genome position (bp)")

ggsave("cdiff_ilds.png", ildsplot, width = 12, height = 8, units = "in")

##### Make rN - rS graphs, 4B-C and S4
# Read in full haplotype csv to calculate pairwise LD
df <- read.delim("C_diff.txt", sep = ",")
max_c <- ncol(df) # index for last haplotype column, as we'll add more later

# Calculate allele frequency, since we only care about common variants
df <- df %>%
  rowwise() %>%
  mutate(allele_freq = mean(c_across(GUT_GENOME142279:GUT_GENOME142432)))
df_common <- subset(df, allele_freq >= 0.2 & site_type != "NC")

mut_types <- c("nonsyn", "syn")

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
# function to compute r2
compute_ld_matrix <- function(geno_df) {
  cor_long <- expand.grid(snp1 = as.numeric(rownames(geno_df)), snp2 = as.numeric(rownames(geno_df)))
  cor_long <- subset(cor_long, snp1 < snp2)
  cor_long$dist <- cor_long$snp2 - cor_long$snp1
  cor_long <- subset(cor_long, dist < 10^5)
  if (nrow(cor_long) > 10^6){ # if there are two many pairs, take a random subset to speed up computation
    cor_long <- cor_long[sample(nrow(cor_long), 10^6), ]
  }
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

linkage_df <- data.frame(snp1 = NA, snp2 = NA, dist = NA, D = NA, r2 = NA, D_prime = NA, site_type = NA)
for (mut in mut_types){
    df_sub <- subset(df_common, site_type == mut)
    if (nrow(df_sub) < 2){
      next
    } else {
      geno_df <- df_sub[,5:max_c]
      rownames(geno_df) <- make.names(df_sub$site_pos, unique = T)
      rownames(geno_df) <- gsub("^X", "", rownames(geno_df))
      
      cor_long <- compute_ld_matrix(geno_df)
      cor_long$site_type <- mut

      linkage_df <- rbind(linkage_df, cor_long)
    }
}

# get rid of any NA values
linkage_df <- na.omit(linkage_df)
linkage_df <- subset(linkage_df, dist < 10^4.5)

# Now create rN - rS graphs for all peak centers
c_diff_sweep_centers <- c(200339, 312063, 611242, 956154, 2740747, 2976139)

mean_r2_df_full <- data.frame(max_dist = NA, mean_dist = NA, num_nonsyn = NA, num_syn = NA, rN = NA, rS = NA, delta_r2 = NA, sweep_center = NA)

for (center in c_diff_sweep_centers){
  
  sub_linkage_df <- subset(linkage_df, snp1 >= (center - dec_dist) & snp2 <= (center + dec_dist))
  
  #dist_bins <- c(0, 6, 40, 100, 500, 1500, 3000, 10000, 32000)
  dist_bins <- c(0, 30, 100, 300, 1000, 10000, 32000)
  
  mean_r2_df <- data.frame(max_dist = dist_bins[2:length(dist_bins)], mean_dist = NA, num_nonsyn = NA, num_syn = NA, rN = NA, rS = NA, delta_r2 = NA)
  for (i in 1:nrow(mean_r2_df)){
    min_dist <- dist_bins[i]
    max_dist <- dist_bins[i+1]
    df_sub <- subset(sub_linkage_df, dist > min_dist & dist <= max_dist)
    
    nonsyn <- subset(df_sub, site_type == "nonsyn")
    syn <- subset(df_sub, site_type == "syn")
    mean_r2_df$mean_dist[i] <- mean(c(min_dist, max_dist))
    mean_r2_df$num_nonsyn[i] <- nrow(nonsyn)
    mean_r2_df$num_syn[i] <- nrow(syn)
    mean_r2_df$rN[i] <- mean(nonsyn$r2, na.rm = T)
    mean_r2_df$rS[i] <- mean(syn$r2, na.rm = T)
    mean_r2_df$delta_r2[i] <- mean(nonsyn$r2, na.rm = T) - mean(syn$r2, na.rm = T)
  }
  
  mean_r2_df$sweep_center <- center
  
  mean_r2_df_full <- rbind(mean_r2_df_full, mean_r2_df)
}

# Supplemental Fig 4 with all sweeps
all_rnrs <- ggplot(data = na.omit(mean_r2_df_full), 
       aes(x = mean_dist, y = delta_r2)) +
  geom_line(color = "blue", linewidth = 2) +
  geom_point(color = "blue", size = 3) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 18)) +
  xlab("Distance between variants (bp)") +
  ylab(expression(rN^2-rS^2)) +
  facet_wrap(~ sweep_center, nrow = 1) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ggsave("cdiff_all_rnrs.pdf", all_rnrs, width = 15, height = 4, units = "in")

# Fig 4B
pilz <- ggplot(data = subset(mean_r2_df_full, sweep_center == 2740747), 
       aes(x = mean_dist, y = delta_r2)) +
  geom_line(color = "blue", linewidth = 2) +
  geom_point(color = "blue", size = 3) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 18)) +
  xlab("Distance between variants (bp)") +
  ylab(expression(rN^2-rS^2)) +
  coord_cartesian(xlim = c(1, 10^4.5),
                  ylim = c(-0.2, 0.2))
ggsave("cdiff_pilz.png", pilz, width = 5, height = 4, units = "in")

# Fig 4C
slayer <- ggplot(data = subset(mean_r2_df_full, sweep_center == 2740747), 
       aes(x = mean_dist, y = delta_r2)) +
  geom_line(color = "blue", linewidth = 2) +
  geom_point(color = "blue", size = 3) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 18)) +
  xlab("Distance between variants (bp)") +
  ylab(expression(rN^2-rS^2)) +
  coord_cartesian(xlim = c(1, 10^4.5),
                  ylim = c(-0.2, 0.2))
ggsave("cdiff_slayer.png", slayer, width = 5, height = 4, units = "in")

##### Make pi scan plot, 4D
# Read in output from pi_sliding_windows.py
cdiff_pi <- read_tsv("C_diff_pi_windows.tsv")

# Merge with significant genes to get colors to match iLDS plot
sig_genes <- ilds[, c("gene_id", "significance")]
sig_genes <- sig_genes %>%
  group_by(gene_id) %>%
  summarise(significance = ifelse(any(significance == "True"), "True", "False"),
            .groups = "drop")
cdiff_pi <- merge(cdiff_pi, sig_genes, by = "gene_id")

# Grey rectangles around two focal sweeps
rects_pi <- data.frame(
  xmin = c(sweep_centers[1]-dec_dist, sweep_centers[2]-dec_dist),
  xmax = c(sweep_centers[1]+dec_dist, sweep_centers[2]+dec_dist),
  ymin = c(0, 0),
  ymax = c(0.165, 0.165),
  fill_color = c("grey", "grey")
)

piplot <- ggplot() + 
  geom_rect(data = rects_pi, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "lightgrey", alpha = 0.5) +
  geom_point(data = subset(cdiff_pi, significance == "False"), 
             aes(x = midpt, y = pi), color = "lightblue") + 
  geom_point(data = subset(cdiff_pi, significance == "True"), 
             aes(x = midpt, y = pi), color = "orange") + 
  theme_bw() + 
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 18)) +
  ylab(expression(Nucleotide~diversity~(pi))) +
  xlab("Genome position (bp)")

ggsave("cdiff_pi.png", piplot, width = 12, height = 8, units = "in")

##### Make haplotype plots, 4E-F
c_diff_sweep_centers <- c(200339, 312063, 611242, 956154, 2740747, 2976139)
haplo_dist <- 10000

# Function to cluster haplotypes according to the distance to the most common haplogroup
get_haplotype_clusters <- function(raw_haps){
  
  # get hamming (manhattan) distances between haplotypes, 
  # and do hierarchical clustering (average ~ UPGMA)
  
  matrix_haps <- t(raw_haps[ ,c(5:139)])
  
  hap_dists <- dist(matrix_haps, method="manhattan")
  hap_clust <- hclust(hap_dists, method="average")
  
  # get haplotype clusters that cluster perfectly (h=0)
  raw_cluster <- cutree(hap_clust, h = 0)
  cluster_counts <- table(raw_cluster) %>% as_tibble() %>%
    arrange(desc(n)) %>%
    mutate(raw_cluster = as.integer(raw_cluster),
           ranked_hap = row_number())
  
  # map samples to clusters, ranked by frequency
  cluster_map <- tibble(sample = names(raw_cluster), raw_cluster) %>%
    left_join(cluster_counts,by='raw_cluster') %>%
    arrange(ranked_hap)
  
  # get clusters' distance to most frequent haplotype cluster
  dist_mat <- as.matrix(hap_dists, labels=TRUE)
  #focal_ben <- rownames(subset(as.data.frame(t(subset(raw_haps, site_pos == 5000))), V1 == 1))
  rep_hap <- cluster_map %>% filter(ranked_hap == 1)
  rep_hap <- rep_hap$sample[1]
  
  #rep_hap <- subset(cluster_map, sample %in% focal_ben)$sample[1]
  
  dist_map <- cluster_map %>% mutate(dist_to_rep = dist_mat[rep_hap,sample]) %>%
    group_by(ranked_hap) %>%
    summarize(ranked_hap=ranked_hap[1],dist_to_rep=dist_to_rep[1]) %>%
    mutate(ranked_hap_dist = rank(dist_to_rep,ties.method="first"))
  
  # join all together
  final_cluster_map <- cluster_map %>% left_join(dist_map, by="ranked_hap") %>%
    rename(ranked_hap_freq = ranked_hap)
  
  return(final_cluster_map)
}

# Pilz domain, Fig 4E
ben_loc <- c_diff_sweep_centers[5]
pilz_df <- subset(df_common, site_pos >= (ben_loc-haplo_dist) & site_pos <= (ben_loc+haplo_dist))

pilz_long <- pivot_longer(pilz_df, cols = GUT_GENOME142279:GUT_GENOME142432,
                           names_to = "sample", values_to = "state")

pilz_clusts <- get_haplotype_clusters(pilz_df)

pilz_long <- merge(pilz_long, pilz_clusts, by = "sample", all = TRUE)
site_pos_all <- sort(unique(pilz_long$site_pos))

haplo_pilz_plot <- ggplot(data = subset(pilz_long)) +
  geom_tile(aes(x = as.factor(site_pos), 
                y = interaction(sample, ranked_hap_dist), 
                fill = interaction(state, site_type))) +
  scale_fill_manual(values = c("1.nonsyn" = "#ed6d52ff", 
                               "1.syn" = "#6f94e6ff", 
                               "0.nonsyn" = "grey", 
                               "0.syn" = "grey"),
                    na.value = "white") +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(breaks = c(site_pos_all[10],
                              site_pos_all[length(site_pos_all)/2],
                              site_pos_all[length(site_pos_all)-10]),
                   labels = c("2.731", "2.734", "2.749")) +
  theme_bw() +
  xlab("Genome position (bp)") +
  ylab("Sample") +
  guides(fill="none") +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none")
ggsave("cdiff_haplo_sweep1.png", haplo_pilz_plot, width = 5, height = 4, units = "in")

# S-layer cassette, Fig 4F
s_ben_loc <- c_diff_sweep_centers[6]
s_df <- subset(df_common, site_pos >= (ben_loc-haplo_dist) & site_pos <= (ben_loc+haplo_dist))

s_long <- pivot_longer(s_df, cols = GUT_GENOME142279:GUT_GENOME142432,
                           names_to = "sample", values_to = "state")

s_clusts <- get_haplotype_clusters(s_df)

s_long <- merge(s_long, s_clusts, by = "sample", all = TRUE)
site_pos_all <- sort(unique(s_long$site_pos))

haplo_s_plot <- ggplot(data = subset(s_long)) +
  geom_tile(aes(x = as.factor(site_pos), 
                y = interaction(sample, ranked_hap_dist), 
                fill = interaction(state, site_type))) +
  scale_fill_manual(values = c("1.nonsyn" = "#ed6d52ff", 
                               "1.syn" = "#6f94e6ff", 
                               "0.nonsyn" = "grey", 
                               "0.syn" = "grey"),
                    na.value = "white") +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(breaks = c(site_pos_all[10],
                              site_pos_all[length(site_pos_all)/2],
                              site_pos_all[length(site_pos_all)-10]),
                   labels = c("2.966", "2.974", "2.985")) +
  theme_bw() +
  xlab("Genome position (bp)") +
  ylab("Sample") +
  guides(fill="none") +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none")
ggsave("cdiff_haplo_sweep2.png", haplo_s_plot, width = 5, height = 4, units = "in")

