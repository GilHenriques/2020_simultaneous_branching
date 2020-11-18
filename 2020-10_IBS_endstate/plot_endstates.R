# Read in 2020-10 branching angles
# Create circular histogram
library(tidyverse)
library(data.table)

`%notin%` <- Negate(`%in%`)

rds_files <- list.files("2020-10_IBS_endstate/", pattern = "rds")

data <- vector(mode = "list", length = length(rds_files))
for (i in seq_along(rds_files)){
  data[[i]] <- readRDS(paste0("2020-10_IBS_endstate/", rds_files[i])) %>% 
  rbindlist(use.names = TRUE) %>% as_tibble()
}

df <- data %>% rbindlist()

# Plot
df <- df %>% mutate(angle_approach = case_when(angle_approach == "90" ~ "pi/2",
                                         angle_approach == "135" ~ "3*pi/4",
                                         angle_approach == "180" ~ "pi",
                                         angle_approach == "225" ~ "-3*pi/4"))

df$angle_approach <-  factor(df$angle_approach, levels = c( "-3*pi/4", "pi/2", "3*pi/4", "pi"))


df %>% 
  separate(class, into = c("Outcome", "outcome2"), sep = " \\(")  %>% 
  ggplot(aes(x = angle_approach, fill = Outcome)) + 
  geom_bar(position = "stack", alpha = 0.9) +
  cowplot::theme_minimal_hgrid() +
  scale_fill_grey() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(labels = parse(text = c( expression(-frac(3*pi, 4)),
                                          expression(frac(pi, 2)),
                                          expression(frac(3*pi, 4)),
                                          expression(pi)
  ))) +
  xlab("Direction of approach") + ylab("Number of replicates") +
  theme(legend.position = "top") 

ggsave("./2020-10_IBS_endstate/endstate_plot.pdf", width = 6, height = 4)
