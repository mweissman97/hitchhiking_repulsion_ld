#### Create figure 5: model of introgression of multilocus adaptation
# Load modules
library(tidyverse)

# Read in file
epdf <- read.csv("epistasis_mean_r2_global.csv")

# Use coarser distance bins for smoother visualization
coarse_dist_bins <- c(10, 30, 150, 500, 1000, 2000, 3000, 4000, 5000, 7500, 10000)
epdf$coarse_dist_bin <- coarse_dist_bins[ findInterval(epdf$dist_bin, coarse_dist_bins, rightmost.closed = TRUE) + 1]
epdf_coarse <- epdf %>%
  group_by(site_type, allele_freq, model, num_sites, coarse_dist_bin) %>%
  summarise(
    mean_r2_weighted = sum(mean_r2 * num_pairs) / sum(num_pairs),
    total_pairs      = sum(num_pairs),
    .groups = "drop"
  )

# Pivot wider to calculate rN - rS
wide_df <- pivot_wider(subset(epdf_coarse, site_type %in% c("nonsyn", "syn")), 
                       id_cols = c(allele_freq, coarse_dist_bin, model, num_sites),
                       names_from = site_type, values_from = mean_r2_weighted)
wide_df$delta_r2 <- wide_df$nonsyn - wide_df$syn

# I modeled several different scenarios, but in the paper we only focus on 3. Take subset of dataframe to streamline visualization
simple_sweep <- subset(wide_df, num_sites == 1 & model == 4)
simple_sweep$model <- "single locus"
simplified_df <- rbind(simple_sweep, 
                    subset(wide_df, num_sites == 10 & model == 4), # just the additive hitchhiking model
                    subset(wide_df, num_sites == 10 & model == 3)) # just the sign epistasis model

# Make figure
epfig <- ggplot(data = subset(simplified_df, 
                     allele_freq == "common")) +
  geom_hline(yintercept = 0) +
  geom_line(aes(x = coarse_dist_bin, y = delta_r2, color = as.factor(model), linetype = as.factor(num_sites)), linewidth = 1) +
  scale_color_manual(values = c("4" = "#dd513a",
                                "3" = "#420a68",
                                "single locus" = "#fca50a"),
                     labels = c("3" = "Epistasis",
                                "4" = "Additive",
                                "single locus" = "Single locus"), 
                     name = "Fitness Model") +
  scale_linetype_manual(values = c("10" = 1, "1" = 4), name = "Number of Loci") +
  theme_bw() +
  scale_x_log10() +
  theme(text = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #legend.position = "none"
  ) +
  xlab("Distance between variants (bp)") +
  ylab(expression(rN^2~-~rS^2)) +
  coord_cartesian(ylim = c(-0.15, 0.15))
ggsave("ep_introgression.png", epfig, width = 6, height = 4, units = "in")


num_site_plot <- ggplot(data = subset(wide_df, 
                     allele_freq == "common" & model != 1 & model != 2)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(x = coarse_dist_bin, 
                y = delta_r2, 
                color = as.factor(num_sites), 
                linetype = as.factor(num_sites)), 
            linewidth = 1) +
  scale_color_manual(values = c("10" = "#420a68", 
                                "2" = "#dd513a",
                                "1" = "#fca50a"), 
                     name = "Number of Loci") +
  scale_linetype_manual(values = c("10" = 1, "2" = 2,"1" = 4), name = "Number of Loci") +
  theme_bw() +
  scale_x_log10() +
  theme(text = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank(), strip.text = element_blank(),
  ) +
  xlab("Distance between variants (bp)") +
  ylab(expression(rN^2~-~rS^2)) +
  coord_cartesian(ylim = c(-0.15, 0.15)) +
  facet_wrap(~model) 
ggsave("num_site_plot.png", num_site_plot, width = 7.5, height = 4, units = "in")

