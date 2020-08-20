# Read in 2020-07 branching angles
# Create circular histogram
library(tidyverse)
library(data.table)

`%notin%` <- Negate(`%in%`)

angles <- readRDS("2020-08 IBS asymmetry/angles_rds")
output <- rbindlist(angles, use.names = TRUE) %>% as_tibble()

plot_angles <- output %>% 
  filter(end_branching) %>% 
  filter(nr_clusters == 2)

plot_angles %>% 
  group_by(alpha) %>% 
  summarize(replicates = n())

plot_angles %>% 
  group_by(alpha) %>% 
  sample_n(600, replace = FALSE) -> plot_angles

angles_2 <- plot_angles %>% mutate(angle = if_else(angle < 0, angle + pi, angle - pi))

angles_origin <- tibble(angle_approach = "225", 
                        origin_angle = -3*pi/4)

pde_angles <- readRDS("2020-07 PDE asymmetry/angle_vs_y_for_comparing_with_ibs_rds")



(plot_ibs <- plot_angles %>% 
  ggplot(aes(x = angle, y = ..count../(sum(..count..)))) + 
  facet_wrap(~ alpha, nrow = 1, labeller = label_bquote(alpha==.(alpha))) +
  geom_histogram(binwidth = diff(seq(from = -pi, to = 0, length = 15))[1],
                 colour = "lightgray", 
                 boundary = pi, 
                 size = 0.1) +
  geom_histogram(data = angles_2, aes(x = angle, y = ..count../(sum(..count..))),
                 binwidth = diff(seq(from = -pi, to = 0, length = 15))[1],
                 fill = "lightgray", 
                 colour = "darkgray",
                 boundary = 0, 
                 size = 0.1) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(-pi, pi),
                     breaks = c(-3*pi/4, -pi/2, -pi/4, 0, pi/4, pi/2, 3*pi/4, pi)
  ) +
  coord_polar(start = pi/2, direction = -1) +
  geom_segment(data = angles_origin,
               aes(x = origin_angle, xend = origin_angle, y = 0.07, yend = 0.01),
               arrow = arrow(length=unit(2, "mm")),
               size = 2.5, color = "white") +
  geom_segment(data = angles_origin,
               aes(x = origin_angle, xend = origin_angle, y = 0.07, yend = 0.01),
               arrow = arrow(length=unit(2, "mm")),
               size = 0.5) +
    geom_vline(data = pde_angles,
               aes(xintercept = angle),
               color = "white", size = 1) +
    geom_vline(data = pde_angles,
               aes(xintercept = angle - pi),
               color = "white", size = 1) +
  geom_vline(data = pde_angles,
               aes(xintercept = angle),
            linetype = "dashed") +
    geom_vline(data = pde_angles,
               aes(xintercept = angle - pi),
               linetype = "dashed") +
  theme_minimal() +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank()) +
  xlab(NULL) + ylab("Frequency") +
  geom_text(data = tibble(
    x = c(-3*pi/4, -pi/2, -pi/4, 0, pi/4, pi/2, 3*pi/4, pi),
    y = 0.13,
    labels = c( expression(-frac(3*pi, 4)),
                expression(-frac(pi, 2)),
                expression(-frac(pi, 4)), 
                "0", 
                expression(frac(pi, 4)), 
                expression(frac(pi, 2)),
                expression(frac(3*pi, 4)),
                expression(pi)
    ) ),
    aes(x = x, y = y, label = labels),
    parse = TRUE, size = 3
  ) +
    scale_y_continuous(breaks = c(0, 0.05, 0.1))
)

ggsave("2020-08 IBS asymmetry/IBS_asymmetry_histogram.pdf", width = 6, height = 2)

saveRDS(plot_ibs, file = "2020-08 IBS asymmetry/plot_ibs_asymmetry_rds")



