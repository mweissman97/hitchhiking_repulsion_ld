library(tidyverse)
library(viridis)

# read in csv
df <- read.csv("paramspace_ld.csv")

# create a wide data frame to calculate rN-rS
df_wide <- pivot_wider(df, 
                       id_cols = c("allele_freq", "dist_bin", "s_ben", "theta", "s_del", "rho_exp", "mu_exp", "h_ben", "h_del", "end_freq", "model", "num_sites"),
                       names_from = "site_type",
                       values_from = "mean_r2")
df_wide$delta_r2 <- df_wide$nonsyn - df_wide$syn

# A) beneficial selection coefficient (s_ben)
sben_plot <- ggplot(data = subset(df_wide, 
                                allele_freq == "common" & 
                                  h_del == h_ben & h_ben == 0.5 &
                                  s_del == 0.001 &
                                  theta == 0.01 &
                                  rho_exp == 6 &
                                  mu_exp == 6 &
                                  end_freq == 0.5 &
                                  model == 0 & num_sites == 1)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(x = dist_bin, y = delta_r2, color = as.factor(s_ben)), linewidth = 1) +
  scale_color_viridis(discrete = T, name = expression(s[ben])) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  ylab(expression(rN^2~-~rS^2)) +
  xlab("Distance between variants (bp)") +
  coord_cartesian(ylim = c(-0.175, 0.175)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 18),
        strip.background = element_blank(), strip.text = element_blank(),
        aspect.ratio = 1); sben_plot
ggsave("plot_sben.png", sben_plot, width = 6, height = 5, units = "in")

# B) deleterious selection coefficient (s_del)
sd_plot <- ggplot(data = subset(df_wide, 
                                allele_freq == "common" & 
                                  h_del == h_ben & h_ben == 0.5 &
                                  s_ben == 0.05 &
                                  theta == 0.01 &
                                  rho_exp == 6 &
                                  mu_exp == 6 &
                                  end_freq == 0.5 &
                                  model == 0 & num_sites == 1)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(x = dist_bin, y = delta_r2, color = as.factor(-s_del)), linewidth = 1) +
  scale_color_viridis(discrete = T, name = expression(s[del])) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  ylab(expression(rN^2~-~rS^2)) +
  xlab("Distance between variants (bp)") +
  coord_cartesian(ylim = c(-0.15, 0.15)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 18),
        strip.background = element_blank(), strip.text = element_blank(),
        aspect.ratio = 1); sd_plot
ggsave("plot_sdel.png", sd_plot, width = 6, height = 5, units = "in")

# C) background mutation rate (mu)
mu_plot <- ggplot(data = subset(df_wide, 
                                 allele_freq == "common" & 
                                   h_del == h_ben & h_ben == 0.5 &
                                   s_ben == 0.05 &
                                   s_del == 0.001 &
                                   theta == 0.01 &
                                   rho_exp == 6 &
                                   end_freq == 0.5 &
                                   model == 0 & num_sites == 1)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(x = dist_bin, y = delta_r2, color = as.factor(-mu_exp)), linewidth = 1) +
  scale_color_viridis(discrete = T, name = expression(mu),
                      labels = parse(text = paste0("10^", levels(factor(-df_wide$mu_exp))))) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  ylab(expression(rN^2~-~rS^2)) +
  xlab("Distance between variants (bp)") +
  coord_cartesian(ylim = c(-0.4, 0.4)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 18),
        strip.background = element_blank(), strip.text = element_blank(),
        aspect.ratio = 1); mu_plot
ggsave("plot_mutrate.png", mu_plot, width = 6, height = 5, units = "in")

# D) recombination rate (rho)
rho_plot <- ggplot(data = subset(df_wide, 
                              allele_freq == "common" & 
                                h_del == h_ben & h_ben == 0.5 &
                                s_ben == 0.05 &
                                s_del == 0.001 &
                                theta == 0.01 &
                                mu_exp == 6 &
                                end_freq == 0.5 &
                                model == 0 & num_sites == 1)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(x = dist_bin, y = delta_r2, color = as.factor(-rho_exp)), linewidth = 1) +
  scale_color_viridis(discrete = T, name = expression(rho),
                      labels = parse(text = paste0("10^", levels(factor(-df_wide$rho_exp))))) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  ylab(expression(rN^2~-~rS^2)) +
  xlab("Distance between variants (bp)") +
  coord_cartesian(ylim = c(-0.15, 0.15)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 18),
        strip.background = element_blank(), strip.text = element_blank(),
        aspect.ratio = 1); rho_plot
ggsave("plot_recomb.png", rho_plot, width = 6, height = 5, units = "in")

# E) dominance coefficient (h)
hplot <- ggplot(data = subset(df_wide, 
                              allele_freq == "common" & 
                                h_del == h_ben &
                                s_ben == 0.05 &
                                s_del == 0.001 &
                                theta == 0.01 &
                                rho_exp == 6 &
                                mu_exp == 6 &
                                end_freq == 0.5 &
                                model == 0 & num_sites == 1)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(x = dist_bin, y = delta_r2, color = as.factor(h_ben)), linewidth = 1) +
  scale_color_viridis(discrete = T, name = "h") +
  theme_bw() +
  xlab("Distance between variants (bp)") +
  ylab(expression(rN^2~-~rS^2)) +
  scale_x_continuous(trans = "log10") +
  coord_cartesian(ylim = c(-0.2, 0.2)) +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 18),
        panel.grid.major = element_blank(),
        strip.background = element_blank(), strip.text = element_blank(),
        aspect.ratio = 1); hplot
ggsave("plot_dominance.png", hplot, width = 6, height = 5, units = "in")

# F) beneficial mutation rate (theta)
type <- ggplot(data = subset(df_wide, 
                             allele_freq == "common" & 
                               h_del == h_ben & h_ben == 0.5 &
                               s_ben == 0.05 &
                               s_del == 0.001 &
                               rho_exp == 6 &
                               mu_exp == 6 &
                               end_freq == 0.5 &
                               model == 0 & num_sites == 1)) +
  geom_line(aes(x = dist_bin, y = delta_r2, color = as.factor(theta)), linewidth = 1) +
  scale_color_viridis(discrete = T, name = expression(theta)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  ylab(expression(rN^2~-~rS^2)) +
  xlab("Distance between variants (bp)") +
  coord_cartesian(ylim = c(-0.15, 0.15)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 18),
        strip.background = element_blank(), strip.text = element_blank(),
        aspect.ratio = 1); type
ggsave("plot_sweeptype.png", type, width = 6, height = 5, units = "in")
