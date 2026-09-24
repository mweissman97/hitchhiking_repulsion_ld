library(tidyverse)
library(viridis)

# part 1: repulsion vs. hitchhiking as a function of distance to beneficial variant
df <- read.csv("params_0.05.0.001.0.01.6.6.0.5.0.5.0.5_all_pairs.csv")
df_common <- subset(df, afreq_bin == "common")

df_common$dist_ben_1 <- abs(df_common$snp1 - 5000)
df_common$dist_ben_2 <- abs(df_common$snp2 - 5000)

df_common$site_type <- ifelse(df_common$snp1 == 5000,
                              "nonsyn_ben",
                              ifelse(df_common$snp2 == 5000,
                                     "nonsyn_ben", df_common$site_type))

df_common$same_side <- ifelse(df_common$snp1 < 5000 & df_common$snp2 < 5000,
                              "same",
                              ifelse(df_common$snp1 > 5000 & df_common$snp2 > 5000,
                                     "same",
                                     "different"))

dist_bins <- c(10, 25, 50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 7500, 10000)
df_common$dist_bin <- dist_bins[ findInterval(df_common$dist, dist_bins, rightmost.closed = TRUE) + 1 ]
df_common_sub <- subset(df_common, same_side == "same" & site_type != "nonsyn_ben")

all_ben_dist_thresh <- c(15, 25, 50, 75)
spat_mean_r2 <- data.frame(dist_ben = NA, site_type = NA, dist_bin = NA, mean_r2 = NA, ben_dist_thresh = NA)
for (thresh in all_ben_dist_thresh){
  df_common_sub$dist_ben <- ifelse(df_common_sub$dist_ben_1 <= thresh & df_common_sub$dist_ben_2 <= thresh,
                                          "near",
                                          ifelse(df_common_sub$dist_ben_1 > thresh & df_common_sub$dist_ben_2 > thresh,
                                                 "far", "mixed"))
  r2_df <- df_common_sub %>%
    group_by(dist_ben, site_type, dist_bin) %>%
    summarize(mean_r2 = mean(r2, na.rm = T),
              ben_dist_thresh = factor(thresh))
  spat_mean_r2 <- rbind(spat_mean_r2, r2_df)
}

spat_mean_r2 <- na.omit(spat_mean_r2)
spat_mean_wide <- pivot_wider(spat_mean_r2,
                              id_cols = c("dist_ben", "dist_bin", "ben_dist_thresh"),
                              names_from = "site_type",
                              values_from = "mean_r2")
spat_mean_wide$delta_r2 <- spat_mean_wide$nonsyn - spat_mean_wide$syn

# figure S8
spat_multi <- ggplot(data = subset(spat_mean_wide, dist_ben != "mixed"), 
       aes(x = dist_bin,
           y = delta_r2,
           color = dist_ben)) +
  geom_hline(yintercept = 0,
             color = "grey") +
  geom_path(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#57106e", "#f98e09"),
                     labels = c("Far", "Near"),
                     name = "") +
  facet_wrap(~ ben_dist_thresh, nrow = 1) +
  scale_x_log10() +
  ylab(expression(rN^2~-~rS^2)) +
  xlab("Distance between variants (bp)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_blank()
        ) 
ggsave("ben_dist_thresh.pdf", spat_multi, width = 10, height = 5, units = "in")

# figure 3A
spat_15 <- ggplot(data = subset(spat_mean_wide, dist_ben != "mixed" & ben_dist_thresh == 15), 
       aes(x = dist_bin,
           y = delta_r2,
           color = dist_ben)) +
  geom_hline(yintercept = 0,
             color = "grey") +
  geom_path(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#57106e", "#f98e09"),
                     labels = c("Far (> 15 bp)", "Near (< 15 bp)"),
                     name = "") +
  scale_x_log10() +
  ylab(expression(rN^2~-~rS^2)) +
  xlab("Distance between variants (bp)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_blank()
        ) 
ggsave("hitch_rep.png", spat_15, width = 8, height = 5, units = "in")

# part 2: repulsion pre vs post sweep
rep_cor_haplos <- read.csv("rep_cor_haplos.csv")
rep_cor_haplos$post_cor_sign <- ifelse(rep_cor_haplos$post_cor < -0.25, "del",
                                       ifelse(rep_cor_haplos$post_cor > 0.25, "pos", "int"))
rep_cor_haplos$pre_cor_sign <- ifelse(rep_cor_haplos$pre_cor < -0.25, "negative correlation",
                                                           ifelse(rep_cor_haplos$pre_cor > 0.25, 
                                                                  "positive correlation",
                                                                  "intermediate correlation"))
rep_cor_haplos$pre_cor_sign[rep_cor_haplos$pre_cor == -1.5] <- "not common"
rep_cor_haplos$pre_cor_sign <- factor(rep_cor_haplos$pre_cor_sign, levels = c("not common", "negative correlation", "intermediate correlation", "positive correlation"))

# figure 3B
rep_hap <- ggplot(data = rep_cor_haplos) +
  geom_histogram(aes(x = post_cor, fill = pre_cor_sign), bins = 10, position = position_fill()) +
  scale_fill_viridis(discrete = T, 
                     option = "magma",
                     name = "Pre-sweep correlation sign", 
                     labels = c("not common" = "Not present", 
                                "negative correlation" = "Negative (< -0.25)", 
                                "intermediate correlation" = "Intermediate (-0.25 - 0.25)", 
                                "positive correlation" = "Positive (> 0.25)")) +
  theme_bw() +
  xlab("Post-sweep correlation (r)") +
  ylab("Proportion of pairs") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 18)); rep_hap

# supplemental figure 9
rep_hap_counts <- ggplot(data = rep_cor_haplos) +
  geom_histogram(aes(x = post_cor, fill = pre_cor_sign), bins = 10) +
  scale_fill_viridis(discrete = T, 
                     option = "magma",
                     name = "Pre-sweep correlation sign", 
                     labels = c("not common" = "Not present", 
                                "negative correlation" = "Negative (< -0.25)", 
                                "intermediate correlation" = "Intermediate (-0.25 - 0.25)", 
                                "positive correlation" = "Positive (> 0.25)")) +
  theme_bw() +
  xlab("Post-sweep correlation (r)") +
  ylab("Proportion of pairs") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 18)); rep_hap_counts
