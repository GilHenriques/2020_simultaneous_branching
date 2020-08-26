# Read in 2020-07 branching angles
# Create circular histogram
library(tidyverse)
library(data.table)

`%notin%` <- Negate(`%in%`)

angles <- readRDS("2020-07 IBS/angles_rds")
output <- rbindlist(angles, use.names = TRUE) %>% as_tibble()

plot_angles <- output %>% 
  filter(end_branching) %>% 
  filter(nr_clusters == 2)

plot_angles %>% 
  group_by(angle_approach) %>% 
  summarize(replicates = n())

plot_angles %>% 
  group_by(angle_approach) %>% 
  sample_n(600, replace = FALSE) -> plot_angles

angles_2 <- plot_angles %>% mutate(angle = if_else(angle < 0, angle + pi, angle - pi))

angles_origin <- tibble(angle_approach = c("45", "90", "135", "180", "225"), 
                        origin_angle = c(pi/4, pi/2, 3*pi/4, pi, -3*pi/4))

(plot_ibs <- plot_angles %>% 
  ggplot(aes(x = angle)) + 
  facet_wrap(~ as.numeric(angle_approach), nrow = 1) +
  geom_histogram(aes(y=..count../(sum(..count..))),
                 binwidth = diff(seq(from = -pi, to = 0, length = 15))[1],
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
               aes(x = origin_angle, xend = origin_angle, y = 0.03, yend = 0.01),
               arrow = arrow(length=unit(2, "mm")),
               size = 2.5, color = "white") +
  geom_segment(data = angles_origin,
               aes(x = origin_angle, xend = origin_angle, y = 0.03, yend = 0.01),
               arrow = arrow(length=unit(2, "mm")),
               size = 0.5) +
  theme_minimal() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_blank()) +
  xlab(NULL) + ylab("Frequency") +
  geom_text(data = tibble(
    x = c(-3*pi/4, -pi/2, -pi/4, 0, pi/4, pi/2, 3*pi/4, pi),
    y = 0.04,
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
  labs(tag = "A")
)

ggsave("./2020-07 IBS/IBS_histogram_600_reps_per_panel.pdf", width = 9, height = 2)

saveRDS(plot_ibs, file = "./2020-07 IBS/plot_ibs_rds")
