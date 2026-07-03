# Load modules
library(tidyverse)
library(viridis)

#### Figure 1: LD curves for s_ben = [0, 0.05] and s_del = [0, -10^-3]
global_ld_df <- read.csv("base_params_mean_r2_global.csv")
# To smooth out figure, use coarser distance bins
new_dist_bins <- c(15, 50, 100, 300, 600, 1000, 1300, 1600, 2000, 2500, 3000, 4000, 5000, 7500, 10000)
global_ld_df$coarse_bin <- new_dist_bins[ findInterval(global_ld_df$dist_bin, new_dist_bins, rightmost.closed = TRUE) + 1 ]

long_df_coarse <- global_ld_df %>%
  group_by(site_type, allele_freq, s_ben, theta, s_del, rho_exp, mu_exp, h_ben, h_del, coarse_bin) %>%
  summarise(
    mean_r2_weighted = sum(mean_r2 * num_pairs) / sum(num_pairs),
    total_pairs      = sum(num_pairs),
    .groups = "drop"
  )

# B-C: separate rN and rS curves
f1b <- ggplot(data = subset(long_df_coarse, allele_freq == "common" & s_ben == 0 & s_del == 0.001),
               aes(x = coarse_bin, 
                   y = mean_r2_weighted,
                   color = site_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("nonsyn" = "red", 
                                  "syn" = "blue"),
                     labels = c("nonsyn" = "Nonsynonymous", 
                                "syn" = "Synonymous"),
                     name = "") +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  ylab(expression(Linkage~disequilibrium~(r^2))) +
  xlab("Distance between variants (bp)") +
  theme(panel.grid = element_blank(),
        text = element_text(size = 18),
        strip.background = element_blank(), strip.text = element_blank()) +
  theme(legend.position = "none")
ggsave("f1b_neutral_rnrs.png", rnrs, width = 6, height = 5, units = "in")

f1c <- ggplot(data = subset(long_df_coarse, allele_freq == "common" & s_ben == 0.05 & s_del == 0.001),
               aes(x = coarse_bin, 
                   y = mean_r2_weighted,
                   color = site_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("nonsyn" = "red", 
                                  "syn" = "blue"),
                     labels = c("nonsyn" = "Nonsynonymous", 
                                "syn" = "Synonymous"),
                     name = "") +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  ylab(expression(Linkage~disequilibrium~(r^2))) +
  xlab("Distance between variants (bp)") +
  theme(panel.grid = element_blank(),
        text = element_text(size = 18),
        strip.background = element_blank(), strip.text = element_blank()) +
  theme(legend.position = "inside",
       legend.position.inside = c(0.8, 0.85))
ggsave("f1c_sweep_rnrs.png", rnrs, width = 6, height = 5, units = "in")

# D-E: rN - rS curves
global_df_wide <- pivot_wider(long_df_coarse, 
                              id_cols = c(allele_freq, coarse_bin, s_ben, theta, s_del, rho_exp, mu_exp, h_ben, h_del), 
                              names_from = site_type, 
                              values_from = mean_r2_weighted,
                              values_fn = mean)
global_df_wide$delta_r2 <- global_df_wide$nonsyn - global_df_wide$syn

f1d <- ggplot(data = subset(global_df_wide, allele_freq == "common" & s_ben == 0 & s_del != 0.01),
               aes(x = coarse_bin, 
                   y = delta_r2,
                   color = as.factor(s_del))) +
  geom_hline(yintercept = 0,
             color = "grey",
             linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#440154", "#5ec962"),
                     labels = c("0", expression(-10^-3)),
                     name = expression(s[del])) +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  ylab(expression(rN^2~-~rS^2)) +
  xlab("Distance between variants (bp)") +
  coord_cartesian(ylim = c(-0.2, 0.2)) +
  theme(panel.grid = element_blank(),
        text = element_text(size = 18),
        strip.background = element_blank(), strip.text = element_blank()) +
  theme(legend.position = "none")
ggsave("f1d_nosweep_delta.png", f1d, width = 6, height = 5, units = "in")

f1e <- ggplot(data = subset(global_df_wide, allele_freq == "common" & s_ben == 0.05 & s_del != 0.01),
               aes(x = coarse_bin, 
                   y = delta_r2,
                   color = as.factor(s_del))) +
  geom_hline(yintercept = 0,
             color = "grey",
             linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#440154", "#5ec962"),
                     labels = c("0", expression(-10^-3)),
                     name = expression(s[del])) +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  ylab(expression(rN^2~-~rS^2)) +
  xlab("Distance between variants (bp)") +
  coord_cartesian(ylim = c(-0.2, 0.2)) +
  theme(panel.grid = element_blank(),
        text = element_text(size = 18),
        strip.background = element_blank(), strip.text = element_blank()) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.2))
ggsave("f1e_sweep_delta.png", f1e, width = 6, height = 5, units = "in")
