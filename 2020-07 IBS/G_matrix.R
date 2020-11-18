library(tidyverse)
library(data.table)

`%notin%` <- Negate(`%in%`)

angles <- readRDS("2020-07 IBS/angles_225_frequent_output_rds")
output <- rbindlist(angles, use.names = TRUE) %>% as_tibble()

exclude_from_distribution_data <- output[output$onset,] %>% 
  filter(t < 6000) %>% 
  pull(file)

distribution_data <- output %>% 
  filter(file %notin% exclude_from_distribution_data) %>% 
  filter(nr_clusters == 2,
         angle_approach == "225") %>% 
  #filter(t %in% c(0, 1500, 3000, 4500))
  filter(t %in% c(0, 200, 400, 600))

distribution_data <- rbind(distribution_data,
                           (distribution_data %>% mutate(angle = angle - pi)),
                           (distribution_data %>% mutate(angle = angle + pi))
)

p1 <- distribution_data %>% 
  ggplot(aes(x = angle, y = t, group = t)) +
  ggridges::geom_density_ridges(alpha = 0.3, scale = 1.5) +
  geom_vline(xintercept = 3*pi/4) +
  coord_cartesian(xlim = c(0,  pi)) + xlab("G matrix orientation (radians)") +
  ylab("Generations") +
  scale_y_continuous(expand = expansion(add = c(0, 20)),
                     breaks = c(0, 200, 400, 600)) +
  scale_x_continuous(expand = c(0,0),
                     breaks = c(0, pi/4, pi/2, 3*pi/4, pi), 
                     labels = c(0, expression(frac(pi,4)), expression(frac(pi,2)), 
                                expression(frac(3*pi,4)), expression(pi))) +
  cowplot::theme_cowplot() +
  theme(axis.text.y = element_text(vjust = 0.5, size = 11),
        axis.text.x = element_text(vjust = 1, size = 11),
        text = element_text(size = 11))



df <- output %>% 
  filter(end_branching | onset,
         angle_approach == "225",
         nr_clusters == 2) 

angle_before_onset <- filter(output, nr_clusters == 2)[which(filter(output, nr_clusters == 2)$onset)-1, ]$angle
angle_branching <- df$angle[df$end_branching]


regr <- readRDS("2020-07 IBS/regression_line")

angles2 <- tibble(angle_before_onset, angle_branching)

angles2 %>% 
    ggplot(aes(x = angle_before_onset, y = angle_branching)) + 
    cowplot::theme_cowplot() +
    theme(legend.position = "none",
          axis.text.y = element_text(vjust = 0.5, size = 11),
          axis.text.x = element_text(vjust = 1, size = 11),
          text = element_text(size = 11)) +
    geom_point(alpha = 0.5) + 
    geom_ribbon(data = regr, aes(x = x, y = y, ymin = ymin, ymax = ymax), alpha = 0.5, fill = "gold") +
    geom_line(data = regr, aes(x = x, y = y), color = "darkorange", size = 1) +
    scale_x_continuous(expand = c(0,0),
                       breaks = c(0, pi/4, pi/2, 3*pi/4, pi), 
                       labels = c(0, expression(frac(pi,4)), expression(frac(pi,2)), 
                                  expression(frac(3*pi,4)), expression(pi))) +
    scale_y_continuous(expand = c(0,0),
                       breaks = c(0, pi/4, pi/2, 3*pi/4, pi), 
                       labels = c(0, expression(frac(pi,4)), expression(frac(pi,2)), 
                                  expression(frac(3*pi,4)), expression(pi))) +
    coord_cartesian(xlim = c(0,pi), ylim= c(0,pi)) +
    xlab("G matrix orientation (radians)\n100 generations prior to branching onset") +
    ylab("Branching direction (radians)") -> p2


p3 <- readRDS("2020-08 IBS bivariate/plot_ibs_bivariate_rds")

final <- cowplot::plot_grid(cowplot::plot_grid((p1 + labs(title = "G-matrix orientation during\napproach to equilibrium") + theme(plot.title = element_text(size = 10, face = "plain"))),
                                      (p2 + labs(title = "G-matrix orientation predicts\nbranching direction") + theme(plot.title = element_text(size = 10, face = "plain"))), 
                                      labels = c("A", "B"), 
                   rel_widths = c(1.7, 2), nrow = 1,
                   align = "h", axis = "b"),
                   (p3+ labs(title = "Branching direction depends on\ndistribution at branching point") + theme(plot.title = element_text(size = 10), plot.title.position = "plot")), 
                   labels = c("", "C"), rel_widths = c(2.5, 1), nrow = 1)


ggsave(filename = "2020-07 IBS/direction_G_matrix.pdf",plot = final, width = 8.5, height = 4)

