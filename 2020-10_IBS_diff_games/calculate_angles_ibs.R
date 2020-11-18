library(tidyverse)
library(data.table)

data_files <- dir(path = "./2020-10_IBS_diff_games/data/180")

bx <- 1.05; cx <- 0.9; dx <- 1.65; essx <-  - (cx - 1)/(2 * (2*bx - cx*dx))
ay <- 1.19; by <- 1.1; cy <-0.9; dy <- 2.14; essy <- 0.2095

pop_size <- 10000
`%notin%` <- Negate(`%in%`)

get_fitness <- function(x, y, xbar, ybar){
   Bx <- (x + xbar) - bx * (x + xbar)^2
   By <- ay * ((y + ybar) - by * (y + ybar)^2)
   Cx <- cx * (x - dx * x^2)
   Cy <- cy * (y - dy * y^2)
   return((Bx-Cx) + (By-Cy))
}

# Use the lines below if you haven't read any files in yet
#angles <- vector(mode = "list", length = length(data_files)) # use this if starting from scratch
#previous_files <- c() # use this if starting from scratch

# Use the lines below if you are reading in angles that have already been calculated
angles <- readRDS("2020-10_IBS_diff_games/angles_180_rds") # read in angles that have already been calculated since before
angles <- c(angles, vector(mode = "list", length = length(data_files)-length(angles))) # extend list so as to include new (not-yet-analysed) data files
previous_files <- angles %>% rbindlist(use.names = TRUE) %>% pull(file)

for(i in 1:length(data_files)){
   file <- data_files[i]
   
   if (file %notin% previous_files) {
      df <- fread(paste0("./2020-10_IBS_diff_games/data/180/", file), col.names = c("t", "x", "y")) %>%   
         mutate(t = t/pop_size) %>% # convert time to generations
         group_by(t) 
      
      branching_end_time <- suppressMessages(df %>% 
                                                summarize(dip_x = diptest::dip.test(x)$statistic[[1]], # Hartigan's dip test for bimodality
                                                          dip_y = diptest::dip.test(x)$statistic[[1]],
                                                          max_dip = max(dip_x, dip_y)) %>% 
                                                mutate(diff = max_dip - 0.15)
      )# first time point when Harrington's dip statistic is > 0.15
      
      if (max(branching_end_time$diff) >= 0) {
         branching_end_time <- branching_end_time %>% 
            filter(diff >= 0) %>% ungroup() %>% 
            filter(t == min(t)) %>% pull(t) %>% unique()
         # test number of branches at end time, based on a sample of 1500 individuals
         db <- suppressMessages(df %>% 
                                   filter(t == branching_end_time) %>% 
                                   select(x, y) %>% 
                                   sample_n(1500) %>% 
                                   # We use DBSCAN: points within epsilon = ess/10 belong to same cluster; minimum points in a cluster is 100
                                   fpc::dbscan(essx/10, MinPts = 100, method = "raw"))
         
         nr_clusters <- db$cluster %>% unique() %>% length()
         
         angles[[i]]$file <- file
         angles[[i]]$nr_clusters <- nr_clusters
         
         if(nr_clusters != 2) {
            angles[[i]]$t <- NA
            angles[[i]]$angle <- NA
            angles[[i]]$onset <- FALSE
            angles[[i]]$end_branching <- FALSE
         } else {
            
            # obtain time for onset of branching
            onset_time <- suppressMessages(df %>% 
                                              filter(t > 0) %>% 
                                              summarize(
                                                 xbar = mean(x), ybar = mean(y), 
                                                 fitness = get_fitness(x, y, xbar, ybar),
                                                 var_fitness = var(fitness)) %>%  # Variance in fitness
                                              ungroup() %>% 
                                              filter(var_fitness == min(var_fitness)) %>% # onset of branching is when variance in fitness is minimum
                                              pull(t) %>% unique()  
            )
            
            # get G matrix orientation up to onset of branching
            G_angle_up_to_onset <- suppressMessages(df %>% 
                                                       filter(t <= branching_end_time) %>% 
                                                       summarize(t = unique(t),
                                                                 varx = var(x), vary = var(y), cov = cov(x,y),
                                                                 evec_x = - (- varx + vary - sqrt(4*cov^2 + varx^2 - 2*varx*vary + vary^2 ))/(2*cov),
                                                                 evec_y = 1,
                                                                 angle = acos(pracma::dot(c(evec_x, evec_y), 
                                                                                          c(1,0))/(pracma::Norm(c(evec_x, evec_y))*pracma::Norm(c(1,0))))) %>% 
                                                       select(t, angle) %>% 
                                                       mutate(onset = ifelse(t == onset_time, TRUE, FALSE)) %>% 
                                                       mutate(end_branching = if_else(t == branching_end_time, TRUE, FALSE))
            )
            
            angles[[i]]$t <- G_angle_up_to_onset$t
            angles[[i]]$angle <- G_angle_up_to_onset$angle
            angles[[i]]$onset <- G_angle_up_to_onset$onset
            angles[[i]]$end_branching <- G_angle_up_to_onset$end_branching
         }
      } else {
         # nothing happens 
         angles[[i]]$file <- file
         angles[[i]]$nr_clusters <- NA
         angles[[i]]$t <- NA
         angles[[i]]$angle <- NA
         angles[[i]]$onset <- NA
         angles[[i]]$end_branching <- NA
      }
      angles[[i]]$angle_approach <- (str_split(file, pattern = "theta")[[1]][2] %>% str_split(pattern = "_"))[[1]][1]
      print(i)
   }
}
beepr::beep("fanfare")


saveRDS(angles, file = "2020-10_IBS_diff_games/angles_180_rds")

