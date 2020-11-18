library(tidyverse)

x <- seq(from = 0, to = 0.238095, length = 100)
y <- seq(from = 0, to = 0.227273, length = 100)

ax <- 1; bx <- 1.05; cx <- 0.9; dx <- 1.65;
ay <- 1.19; by <- 1.1; cy <- 0.9; dy <- 2.14;

Bx <- ax*(x - bx*(x)^2)
Cx <- cx*(x - dx*x^2)

By <- ay*(x - by*(x)^2)
Cy <- cy*(x - dy*x^2)

df <- bind_rows(
tibble(x = x, y = Bx, Function = "Benefit", Game = "1"),
tibble(x = y, y = By, Function = "Benefit", Game = "2"),
tibble(x = x, y = Cx, Function = "Cost", Game = "1"),
tibble(x = y, y = Cy, Function = "Cost", Game = "2"),
)

payoff_plot <- ggplot(data = df) +
  geom_line(aes(x = x, y = y, color = Game, linetype = Function), size = 1) +
  geom_vline(aes(xintercept = 0.0813)) +
  geom_vline(aes(xintercept = 0.209538), color = "darkgrey") +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  xlab(expression(paste("Investment strategy (", italic(z), ")"))) +
  ylab("Benefit or cost") +
  scale_color_grey() +
  theme(legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 9),
        legend.spacing.y = unit(0.025, "cm"),
        legend.background = element_rect(fill = "transparent", colour = NA), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", colour = NA) # get rid of legend panel bg
        )

saveRDS(payoff_plot, file = "./2020-10_IBS_diff_games/plot_payoff_rds")

payoff_plot <- read_rds("./2020-10_IBS_diff_games/plot_payoff_rds") + 
  labs(subtitle = "Cost and benefit for games 1 and 2") +
  theme(plot.subtitle = element_text(size = 10))

ibs_plot <- read_rds("./2020-10_IBS_diff_games/plot_ibs_rds") + 
  labs(subtitle = "Branching direction with games 1 and 2")+
  theme(plot.subtitle = element_text(size = 10))

random_plot <- read_rds("./2020-11_IBS_random_diff_games/plot_ibs_rds") + 
  labs(subtitle = "Branching direction with game 1 and random games")+
  theme(plot.subtitle = element_text(size = 10))

full_plot <- cowplot::plot_grid(payoff_plot, ibs_plot, random_plot, ncol = 1,
                                rel_heights = c(0.9, 1, 1), labels = c("A", "B", "C"))

ggsave("./2020-10_IBS_diff_games/diff_games_plot.pdf", plot = full_plot, width = 4, height = 6)
   
