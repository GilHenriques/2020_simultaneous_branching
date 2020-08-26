library(tidyverse)
library(data.table)
`%notin%` <- Negate(`%in%`)


angles <- readRDS("2020-08 IBS asymmetry/angles_rds")
output <- rbindlist(angles, use.names = TRUE) %>% as_tibble() %>% filter(nr_clusters == 2)

files_to_exclude <- output %>% 
  filter(onset) %>% 
  filter(t < 4000) %>% 
  pull(file) %>% unique()

df <- output %>% filter(file %notin% files_to_exclude)

rows_onset <- which(df$onset)
rows_1000_before <- which(df$onset)-1
rows_2000_before <- which(df$onset)-2
rows_3000_before <- which(df$onset)-3

df_G_plot <- rbind(
  (df[rows_1000_before, ] %>% select(angle, alpha) %>% mutate(time_before_onset = "1000 generations\n before onset")),
  (df[rows_2000_before, ] %>% select(angle, alpha) %>% mutate(time_before_onset = "2000 generations\n before onset")),
  (df[rows_3000_before, ] %>% select(angle, alpha) %>% mutate(time_before_onset = "3000 generations\n before onset")),
  (df[df$t==0, ] %>% select(angle, alpha) %>% mutate(time_before_onset = "Initial condition"))
  
)

df_G_plot <- bind_rows(df_G_plot,
                       mutate(df_G_plot, angle = angle-pi),
                       mutate(df_G_plot, angle = angle+pi))


G_plot <- df_G_plot %>%   
  ggplot() +
  ggridges::geom_density_ridges(aes(x = angle, y = reorder(time_before_onset, desc(time_before_onset)), fill = as.character(alpha)),
                                alpha = 0.3, scale = 1.5) +
  geom_vline(xintercept = 3*pi/4) +
  scale_fill_viridis_d(option = "C", name = expression(alpha)) +
  cowplot::theme_cowplot() +
  theme(axis.text.y = element_text(vjust = 0), legend.position = "bottom", axis.title.y = element_blank()) +
  scale_y_discrete(expand = expansion(mult = c(0, 0.25))) +
  scale_x_continuous(expand = c(0,0),
                    breaks = c(0, pi/4, pi/2, 3*pi/4, pi), 
                    labels = c(0, expression(frac(pi,4)), expression(frac(pi,2)), 
                               expression(frac(3*pi,4)), expression(pi))) +
  coord_cartesian(xlim = c(0,pi)) + xlab("G matrix orientation (radians)")

saveRDS(G_plot, file = "2020-08 IBS asymmetry/G_plot_rds")

