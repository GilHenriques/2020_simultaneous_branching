# Read in 2020-08 branching angles
# Create circular histogram

library(tidyverse)

df <- readRDS("./2020-08 OSS asymmetry/branching_angles_rds")

angles_origin <- tibble(angle_of_approach = "225", 
                        origin_angle = -3*pi/4)

df_plot <- df %>% 
  drop_na() %>% 
  group_by(asymmetry) %>% 
  sample_n(7950, replace = FALSE)

# Copy of angles for white histogram bars
df_2 <- df_plot %>% mutate(angle_of_branching_2 = if_else(angle_of_branching < 0, angle_of_branching + pi, angle_of_branching - pi), asymmetry2 = as.numeric(asymmetry))

(plot_oss <- df_plot %>% 
    mutate(asymmetry2 = as.numeric(asymmetry)) %>% 
  ggplot(aes(x = angle_of_branching, y = ..count../(sum(..count..)))) + 
  geom_histogram(binwidth = diff(seq(from = 0, to = pi, length = 15))[1], 
                 colour = "lightgray", 
                 boundary = pi, 
                 size = 0.1) +
  geom_histogram(data = df_2, aes(x = angle_of_branching_2, y = ..count../(sum(..count..))),
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
  facet_wrap(.~asymmetry2, nrow = 1, labeller = label_bquote(alpha==.(asymmetry2))) +
  geom_segment(data = angles_origin,
               aes(x = origin_angle, xend = origin_angle, y = 0.04, yend = 0.01),
               arrow = arrow(length=unit(2, "mm")),
               size = 2.5, color = "white") +
  geom_segment(data = angles_origin,
               aes(x = origin_angle, xend = origin_angle, y = 0.04, yend = 0.01),
               arrow = arrow(length=unit(2, "mm")),
               size = 0.5) +
  theme_minimal() +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank()) +
  xlab(NULL) + ylab("Frequency") +
  geom_text(data = tibble(
    x = c(-3*pi/4, -pi/2, -pi/4, 0, pi/4, pi/2, 3*pi/4, pi),
    y = 0.06,
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
    parse = TRUE, size = 3) +
    scale_y_continuous(breaks =c(0, 0.02, 0.04, 0.06))
  ) 
  
ggsave("./2020-08 OSS asymmetry/OSS_histogram_7950_reps_per_panel.pdf", width = 6, height = 2)

saveRDS(plot_oss, file = "./2020-08 OSS asymmetry/plot_oss_asymmetry_rds")
