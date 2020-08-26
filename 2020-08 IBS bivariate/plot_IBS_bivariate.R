library(tidyverse)
library(data.table)

`%notin%` <- Negate(`%in%`)

angles <- readRDS("2020-08 IBS bivariate/angles_rds")
output <- rbindlist(angles, use.names = TRUE) %>% as_tibble()
output <- output %>% mutate(rho = ifelse(alpha == "0.5", 0.5, alpha),
                  rho = ifelse(rho == "Minus0.5", -0.5, rho),
                  rho = as.numeric(rho)) 

plot_angles <- output %>% 
  filter(end_branching) %>% 
  filter(nr_clusters == 2)

plot_angles %>% 
  group_by(rho) %>% 
  summarize(replicates = n())

plot_angles %>% 
  group_by(rho) %>% 
  sample_n(600, replace = FALSE) -> plot_angles

angles_2 <- plot_angles %>% mutate(angle = if_else(angle < 0, angle + pi, angle - pi))

(plot_ibs <- plot_angles %>% 
    ggplot(aes(x = angle, y = ..count../(sum(..count..)))) + 
    facet_wrap(~ rho, nrow = 2, labeller = label_bquote(rho==.(rho))) +
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
    theme_minimal() +
    theme(strip.background = element_blank(),
          axis.text.x = element_blank()) +
    xlab(NULL) + ylab("Frequency") +
    geom_text(data = tibble(
      x = c(-3*pi/4, -pi/2, -pi/4, 0, pi/4, pi/2, 3*pi/4, pi),
      y = 0.1,
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
    ) 
)


ggsave("2020-08 IBS bivariate/IBS_bivariate_histogram.pdf", width = 2.5, height = 4.5)

saveRDS(plot_ibs, file = "2020-08 IBS bivariate/plot_ibs_bivariate_rds")
