library(tidyverse)

covariance <- function(axis, f, threshold){
  f[f < threshold] <- 0
  mean_x <- (axis * rowSums(f)) %>% sum 
  mean_y <- (axis * colSums(f)) %>% sum 
  cov <- 0
  for(index_x in seq_along(axis)){ 
    for(index_y in seq_along(axis)){
      cov <- cov + (axis[index_x] - mean_x)*(axis[index_y] - mean_y)*f[index_x, index_y]
    }
  }
  return(cov)
}

var_x <- function(axis, f, threshold){
  f[f < threshold] <- 0
  mean_x <- (axis * rowSums(f)) %>% sum
  varx <- 0
  for(index_x in seq_along(axis)){
    varx <- varx + (axis[index_x] - mean_x)^2 * rowSums(f)[index_x]
  }
  return(varx)
}

var_y <- function(axis, f, threshold){
  f[f < threshold] <- 0
  mean_y <- (axis * colSums(f)) %>% sum
  vary <- 0
  for(index_y in seq_along(axis)){
    vary <- vary + (axis[index_y] - mean_y)^2 * colSums(f)[index_y]
  }
  return(vary)
}

multimodal <- function(axis, f){ # returns true if the marginal distribution is multimodal along either x or y
  marginal_y <- colSums(f); marginal_x <- rowSums(f)
  multimodal_y <- length(rle(marginal_y > (max(marginal_y))/10)$values) > 3
  multimodal_x <- length(rle(marginal_x > (max(marginal_x))/10)$values) > 3
  return(any(multimodal_x, multimodal_y))
}

paths <- dir(recursive = T, pattern = "full_population")
numbers <- as.numeric(map_chr(paths, function(x) strsplit(strsplit(x, "/")[[1]][2], "_")[[1]][1]))
selection_strength <- as.numeric(map_chr(paths, function(x) strsplit(strsplit(x, "_s=")[[1]][2], "_mu")[[1]][1]))
mu0 <- as.numeric(map_chr(paths, function(x) strsplit(strsplit(x, "_mu=")[[1]][2], "_angle")[[1]][1]))
zeta0 <- as.numeric(map_chr(paths, function(x) strsplit(strsplit(x, "zeta=")[[1]][2], "_s")[[1]][1]))
theta <- as.numeric(map_chr(paths, function(x) strsplit(strsplit(x, "angle=")[[1]][2], "/")[[1]][1]))

df <- tibble(path = paths, numbers, selection_strength, mu0, zeta0, theta) 
threshold <- 1e-5

df$covariance <- NA
df$var_x <- NA
df$var_y <- NA
df$time <- NA
df$multimodal <- NA
pb <- txtProgressBar(min = 0, max = nrow(df), style = 3)
for(i in 1:nrow(df)){
  path <- df[i, ]$path
  f <- readRDS(path)[[1]]
  axis <- readRDS(path)[[2]]
  time <- readRDS(path)[[4]]
  df[i, ]$covariance <- covariance(axis, f, threshold)
  df[i, ]$var_x <- var_x(axis, f, threshold)
  df[i, ]$var_y <- var_y(axis, f, threshold)
  df[i, ]$time <- time
  df[i, ]$multimodal <- multimodal(axis, f)
  setTxtProgressBar(pb, i)
}

df_not_branched <- df %>%
  arrange(zeta0, mu0, selection_strength, time) %>% 
  group_by(zeta0, mu0, selection_strength) %>% 
  mutate(branched = cumsum(multimodal)) %>%
  mutate(branched = if_else(branched == 0, FALSE, TRUE)) %>% 
  filter(branched == FALSE)

saveRDS(df_not_branched, file = "data_not_branched")

df_not_branched %>%  
  filter(theta == 225) %>% # or 135
  ggplot(aes(x = time, y = covariance, color = as.factor(mu0))) +
  geom_line() +
  geom_hline(yintercept = 0) +
  facet_grid(zeta0 ~ selection_strength, scales = "free") +
  theme(legend.position = "bottom")

df_not_branched %>%  
  filter(theta == 270) %>% # or 270
  ggplot(aes(x = time, y = var_x - var_y, color = as.factor(mu0))) +
  geom_line() +
  geom_hline(yintercept = 0) +
  facet_grid(zeta0 ~ selection_strength, scales = "free") +
  theme(legend.position = "bottom")

df_not_branched %>%  
  filter(theta == 225, mu0 == 0.01, zeta0 == 2e-07, selection_strength == 100) %>% # or 135
  ggplot(aes(x = time, y = covariance)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "grey") +
  cowplot::theme_half_open()


# We want smaller mu, smaller selection strength, and and bigger zeta. E.g. mu = 0.001, z = 2e-7, mu = 0.001


