# This script combines the IBS and OSS histograms

library(tidyverse)

plot_oss <- readRDS("./2020-01 OSS/plot_oss_rds") + 
  labs(title = "Branching direction with OSS") +
  theme(plot.title = element_text(size = 10))
plot_ibs <- readRDS("./2020-07 IBS/plot_ibs_rds") + 
  labs(title = "Branching direction with IBS") +
  theme(plot.title = element_text(size = 10))

cowplot::plot_grid((plot_ibs + labs(tag = "")), 
                   (plot_oss + labs(tag = "")), 
                   nrow = 2, labels = c("A", "B"))

ggsave("./2020-07 IBS/OSS_IBS_histogram.pdf", width = 9, height = 4.25)
