#### Distance / frequency matching
neutral_paired_df <- read.csv("nosweep_paired_delta_r2.csv")
sweep_paired_df <- read.csv("sweep_paired_delta_r2.csv")
freq_paired_df <- rbind(neutral_paired_df, sweep_paired_df)
common_freq_df <- subset(freq_paired_df, afreq_bin == "common")
common_freq_df$dist_bin <- cut(common_freq_df$dist, 
                               breaks = c(0, 10, 30, 100, 300, 1000, 1300, 2000, 3000, 4000, 5000, 10000),
                               labels = c(5, 20, 65, 200, 650, 1200, 1650, 2500, 3500, 4500, 7500))

common_freq_mean <- common_freq_df %>% 
  group_by(s_ben, dist_bin) %>%
  summarize(mean_delta_r2 = mean(delta_r2))

common_freq_mean$dist_bin <- as.numeric(as.character(common_freq_mean$dist_bin))

# Neutral plot
freq_paired_neut <- ggplot() +
  geom_point(data = subset(common_freq_df, s_ben == 0), 
             aes(x = dist, y = delta_r2), alpha = 0.05, color = "grey") +
  geom_line(data = subset(common_freq_mean, s_ben == 0),
            aes(x = dist_bin, y = mean_delta_r2), color = "darkgreen", linewidth = 1.5) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  coord_cartesian(ylim = c(-0.25, 0.25), 
                  xlim = c(1, 10000)) +
  xlab("Distance between variants (bp)") +
  ylab(expression(rN^2~-~rS^2)) +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 18),
        strip.background = element_blank(), strip.text = element_blank(),
        legend.position = "none")
ggsave("freqpaired_neut.png", freq_paired_neut, width = 5, height = 4, units = "in")

# Sweep plot
freq_paired_sweep <- ggplot() +
  geom_point(data = subset(common_freq_df, s_ben == 0.05), 
             aes(x = dist, y = delta_r2), alpha = 0.05, color = "grey") +
  geom_line(data = subset(common_freq_mean, s_ben == 0.05),
            aes(x = dist_bin, y = mean_delta_r2), color = "darkgreen", linewidth = 1.5) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  coord_cartesian(ylim = c(-0.25, 0.25), 
                  xlim = c(1, 10000)) +
  xlab("Distance between variants (bp)") +
  ylab(expression(rN^2~-~rS^2)) +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 18),
        strip.background = element_blank(), strip.text = element_blank(),
        legend.position = "none")
ggsave("freqpaired_sweep.png", freq_paired_sweep, width = 5, height = 4, units = "in")
