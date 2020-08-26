library(tidyverse)
library(data.table)

b <- 1.05; c <- 0.9; d <- 1.65; ess <-  - (c - 1)/(2 * (2*b - c*d))
pop_size <- 10000
`%notin%` <- Negate(`%in%`)

get_fitness <- function(x, y, xbar, ybar){
  bx <- (x + xbar) - b * (x + xbar)^2
  by <- (y + ybar) - b * (y + ybar)^2
  cx <- c * (x - d * x^2)
  cy <- c * (y - d * y^2)
  return((bx-cx) + (by-cy))
}

files_ibs <- tibble(file = dir(path = "./2020-07 IBS/data")) %>% 
  mutate(tmp = file) %>% 
  separate(tmp, sep = "theta", into = c("tmp", "angle")) %>% 
  separate(angle, sep = "_", into = c("angle", "tmp2")) %>% 
  select(file, angle)

files_ibs %>% 
  slice(5) %>% 
  pull(file) -> example_ibs

df_ibs <- read_csv(paste0("./2020-07 IBS/data/", example_ibs), col_names = c("t", "x", "y")) %>% 
  mutate(t = t/pop_size) %>% # convert time to generations
  group_by(t) 

df_ibs %>% sample_n(500) %>% ggplot(aes(x=x, y=y)) + geom_point()

ibs_info <- readRDS("./2020-07 IBS/angles_rds") %>% rbindlist(use.names = TRUE) %>% as_tibble() %>% filter(file == example_ibs)

branching_end_time_ibs <- ibs_info %>% filter(end_branching) %>% pull(t)
onset_time_ibs <- ibs_info %>% filter(onset) %>% pull(t)
branching_angle_ibs <- ibs_info %>% filter(end_branching) %>% pull(angle)
onset_angle_ibs <- ibs_info %>% filter(onset) %>% pull(angle)

ibs_at_branching <- df_ibs %>%
  filter(t == branching_end_time_ibs) %>% 
  ungroup() 
db <- ibs_at_branching %>% 
  select(x, y) %>% 
  fpc::dbscan(ess/10, MinPts = 100, method = "raw")
ibs_at_branching$cluster <- db$cluster
branching_line_ibs <- ibs_at_branching %>% 
  group_by(cluster) %>% 
  summarize(x = mean(x), y = mean(y)) 
slope_ibs <- (branching_line_ibs$y[2] - branching_line_ibs$y[1]) / (branching_line_ibs$x[2] - branching_line_ibs$x[1])
intercept_ibs <- branching_line_ibs$y[1] - slope_ibs*branching_line_ibs$x[1]
branching_line_ibs <- tibble(slope = slope_ibs, intercept = intercept_ibs)

onset_line_ibs <- df_ibs %>% 
  filter(t == onset_time_ibs) %>% 
  summarize(mean_x = mean(x), mean_y = mean(y)) %>% 
  mutate(slope = tan(onset_angle_ibs))

(panel_A <- df_ibs %>%
  group_by(t) %>% 
  sample_n(200) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_point(alpha = 0.1) +
  #geom_abline(data = onset_line_ibs, aes(slope = slope, intercept = (mean_y - slope*mean_x)), color = "firebrick") +
  geom_abline(data = branching_line_ibs, aes(slope = slope, intercept = intercept), color = "gold") +
  geom_point(data = (df_ibs %>% filter(t == branching_end_time_ibs)), color = "gold", alpha = 0.03) +
  geom_point(data = (df_ibs %>% filter(t == onset_time_ibs)), color = "firebrick", alpha = 0.03) +
  coord_fixed(xlim = c(0.05, 0.11), ylim = c(0.05, 0.11)) +
  cowplot::theme_cowplot() +
  labs(x = expression(z[1]), y = expression(z[2]))
)

nr_poly <- 10 # nr of invasions with polymorphic population until branching is "complete"

data_oss <- paste0("./2020-01 OSS/data/theta_45_output_list") %>% 
  read_rds() %>% 
  map(. %>% drop_na())

df_oss <- data_oss[[29]] 

df_oss %>%
  ggplot(aes(x = x, y = y)) +
  geom_point() +
  coord_fixed()

time_strains <- df_oss %>% 
  drop_na() %>% 
  group_by(t) %>% 
  summarize(nr_strains = n()) %>%  
  mutate(poly = nr_strains >= 2) 

counter <- 0; time_strains$sequential_poly <- 0;
for(j in 1:nrow(time_strains)){
  if(time_strains[j, ]$poly) { counter <- counter + 1 } else { counter <- 0 }
  time_strains[j, ]$sequential_poly <- counter
}

time_branching <- time_strains %>% 
  filter(sequential_poly == nr_poly) %>% 
  slice(1) %>% pull(t)

branching_line_oss <- df_oss %>% 
  filter(t == time_branching) 
slope_oss <- (branching_line_oss$y[2] - branching_line_oss$y[1]) / (branching_line_oss$x[2] - branching_line_oss$x[1])
intercept_oss <- branching_line_oss$y[1] - slope_oss*branching_line_oss$x[1]
branching_line_oss <- tibble(slope = slope_oss, intercept = intercept_oss)

(oss_plot <- df_oss %>%
  ggplot(aes(x = x, y = y)) +
  geom_point() +
  geom_point(data = filter(df_oss, t == time_branching), color = "gold", size = 3) +
  geom_abline(data = branching_line_oss, aes(slope = slope, intercept = intercept), color = "gold") +
  coord_fixed(ylim = c(0.058, 0.11), xlim = c(0.058, 0.11)) +
  cowplot::theme_cowplot() +
  scale_x_continuous(breaks =  c(0.06, 0.08, 0.10)) +
  scale_y_continuous(breaks =  c(0.06, 0.08, 0.10)) +
  labs(x = expression(z[1]), y = expression(z[2]))
)

ibs_oss_plot <- cowplot::plot_grid(panel_A, oss_plot, nrow = 1, labels = c("A", "B"))

ggsave(filename = "2020-07 IBS/example_replicate.pdf", plot = ibs_oss_plot, width = 5.75, height = 2.75)

