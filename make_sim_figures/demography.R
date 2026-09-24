library(tidyverse)
library(viridis)

simple_comp_df <- read.csv("demographic_contractions_deltar2.csv")

demo <- ggplot(data = subset(simple_comp_df, s_focal == 0)) +
  geom_line(aes(x = dist_bin, y = delta_r2, color = demography, linetype = demography), linewidth = 1) +
  scale_color_viridis(discrete = T, labels = c("long" = "Long, shallow contraction", "short" = "Short, sharp contraction", "none" = "Constant N"),
                      name = "Demography") +
  scale_linetype_manual(values = c("long" = "longdash", "short" = "dotdash", "none" = "solid"),
                        labels = c("long" = "Long, shallow contraction", "short" = "Short, sharp contraction", "none" = "Constant N"),
                        name = "Demography") +
  geom_hline(yintercept = 0) +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  ylab("rN - rS") +
  xlab("Distance between variants (bp)") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 18)) +
  coord_cartesian(ylim = c(-0.16, 0.1)); demo
ggsave("hitch_demography.png", demo, width = 8, height = 5, units = "in")
