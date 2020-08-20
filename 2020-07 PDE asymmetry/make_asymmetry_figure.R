library(tidyverse)
ess <- 0.0813
directories <- dir("2020-07 PDE asymmetry/data_for_comparing_with_ibs/")

frequencies <- tibble(x = double(), y = double(), Freq = double(), snapshot = double(), asymmetry = character())
for(directory in directories){
  asymmetry <- (directory %>% str_split(pattern = "alpha="))[[1]][2] %>% 
    as.numeric() %>% round(3) %>% as.character()
  
  all_paths <- dir(path = paste0("2020-07 PDE asymmetry/data_for_comparing_with_ibs/", directory), pattern = "full_population") %>% 
    gtools::mixedsort()
  #if(as.numeric(asymmetry) < 3) {
  #  paths <- all_paths[seq(from = 1, to = 500/5, length = 5)]
  #} else { 
  #  paths <- all_paths[seq(from = 1, to = 400/5, length = 5)]
  #}
  paths <- all_paths[seq(from = 1, to = 400/3, length = 5)]
  
  axis <- readRDS(paste0("2020-07 PDE asymmetry/data_for_comparing_with_ibs/", directory, "/", paths[1]))[[2]]
  axisY <- axis[axis > 0.04 & axis < 0.115]
  axisX <- axis[axis > 0.04 & axis < 0.115]
  
  for(snapshot in 1:length(paths)) {
    distr <- readRDS(paste0("2020-07 PDE asymmetry/data_for_comparing_with_ibs/", directory, "/", paths[snapshot]))[[1]]
    distr <- distr[axis > 0.04 & axis < 0.115, axis > 0.04 & axis < 0.115]
    colnames(distr) = axisY
    rownames(distr) = axisX
    distr <- as.data.frame.table(distr) %>% 
      tibble() %>% 
      mutate(x = as.numeric(as.character(Var1)), y = as.numeric(as.character(Var2)), snapshot = snapshot, asymmetry = asymmetry) %>% 
      dplyr::select(x, y, Freq, snapshot, asymmetry)
    frequencies <- rbind(distr, frequencies)
  }
  
  print(asymmetry)
}

(
plot_freqs <- frequencies %>%
  ggplot(aes(x = x, y = y, z = Freq)) +
  geom_hline(yintercept = ess) +
  geom_vline(xintercept = ess) +
  geom_contour(data = filter(frequencies, snapshot == 1), alpha = 0.2, color = "black") +
  geom_contour(data = filter(frequencies, snapshot == 2), alpha = 0.4, color = "black") +
  geom_contour(data = filter(frequencies, snapshot == 3), alpha = 0.6, color = "black") +
  geom_contour(data = filter(frequencies, snapshot == 4), alpha = 0.8, color = "black") +
  geom_contour(data = filter(frequencies, snapshot == 5), alpha = 1, color = "black") +
  facet_wrap(.~asymmetry, labeller = label_bquote(alpha==.(asymmetry)),nrow = 1) +
  cowplot::theme_cowplot() +
  cowplot::panel_border() +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        aspect.ratio = 1) +
  labs(x = expression(z[1]), y = expression(z[2]))
)


angles <- readRDS("2020-07 PDE asymmetry/angle_vs_y_rds")

(plot_angles <- angles %>% 
  ggplot(aes(x = as.numeric(alpha), y = angle)) + geom_point() + 
  scale_y_continuous(limits = c(pi/2, pi),
                     breaks = c(pi/2, 2*pi/3, 5*pi/6, pi), 
                     labels = c( expression(frac(pi, 2)),
                                 expression(frac(2*pi, 3)),
                                 expression(frac(5*pi, 6)),
                                 expression(pi)) ) +
    xlab(expression("Asymmetry between traits ("*alpha*")")) +
    ylab("Angle of branching")  +
  cowplot::theme_cowplot()
)





plot_oss <- readRDS("2020-08 OSS asymmetry/plot_oss_asymmetry_rds")
plot_ibs <- readRDS("2020-08 IBS asymmetry/plot_ibs_asymmetry_rds")
plot_G <- readRDS("2020-08 IBS asymmetry/G_plot_rds")

plot_final <- cowplot::plot_grid(
  cowplot::plot_grid(
    cowplot::plot_grid(plot_oss, plot_ibs, nrow = 2, labels = c("A", "B")),
    (plot_G + theme(text = element_text(size = 10), axis.text = element_text(size = 10))), 
    nrow = 1, rel_widths = c(1, 0.6), labels = c(NA, "C")
  ),
  cowplot::plot_grid((plot_angles + theme(text = element_text(size = 10), axis.text = element_text(size = 10))), 
                     (plot_freqs + theme(text = element_text(size = 10), axis.text = element_text(size = 10))), 
                     rel_widths = c(1.5, 2.5), labels = c("D", "E"),
                     align = "h", axis = "bt"),
  ncol = 1, rel_heights = c(1, 0.55)
)

ggsave(filename = "2020-07 PDE asymmetry/figure_asymmetry.pdf", plot = plot_final, width = 9, height = 6)
