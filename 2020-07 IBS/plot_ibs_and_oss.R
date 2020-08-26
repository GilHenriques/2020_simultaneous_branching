# This script combines the IBS and OSS histograms

library(tidyverse)

plot_oss <- readRDS("./2020-01 OSS/plot_oss_rds")
plot_ibs <- readRDS("./2020-07 IBS/plot_ibs_rds")

cowplot::plot_grid((plot_ibs + labs(tag = "")), 
                   (plot_oss + labs(tag = "")), 
                   nrow = 2, labels = c("A", "B"))

ggsave("./2020-07 IBS/OSS_IBS_histogram.pdf", width = 9, height = 4)
