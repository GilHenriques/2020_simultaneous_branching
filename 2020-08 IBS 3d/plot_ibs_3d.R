library(tidyverse)
library(data.table)

`%notin%` <- Negate(`%in%`)

angles <- readRDS("2020-08 IBS 3d/angles_rds")
output <- rbindlist(angles, use.names = TRUE) %>% as_tibble()


plot_angles <- output %>% 
  drop_na() %>% 
  filter(nr_clusters == 2)
# 597 replicates

plot_angles <- plot_angles %>% 
  pivot_longer(-c(file, nr_clusters, t)) %>% 
  mutate(name = ifelse(name == "alpha", "varphi", "psi"))

plot_angles$name <- factor(plot_angles$name, levels = c("varphi", "psi"))


angles_2 <- plot_angles %>% mutate(value = if_else(value < 0, value + pi, value - pi))


(plot_ibs <- plot_angles %>% 
    ggplot(aes(x = value)) + 
    facet_wrap(~ name, nrow = 1, labeller = label_parsed) +
    geom_histogram(aes(y = ..count../(sum(..count..))),
                   binwidth = diff(seq(from = -pi, to = 0, length = 15))[1],
                   colour = "lightgray", 
                   boundary = pi, 
                   size = 0.1) +
    geom_histogram(data = angles_2, aes(x = value, y = ..count../(sum(..count..))),
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
    theme_minimal() +
    theme(strip.background = element_blank(),
          axis.text.x = element_blank()) +
    xlab(NULL) + ylab("Frequency") +
    geom_text(data = tibble(
      x = c(-3*pi/4, -pi/2, -pi/4, 0, pi/4, pi/2, 3*pi/4, pi),
      y = 0.07,
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
    geom_segment(data = filter(plot_angles, name == "varphi"), x = pi/2, xend = pi/2, 
                 y = 0.055, yend = 0.02,
                 arrow = arrow(length=unit(2, "mm")),
                 size = 0.3) 
)
ggsave("2020-08 IBS 3d/IBS_3d.pdf", width = 4, height = 2)


