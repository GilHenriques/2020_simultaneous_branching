library(tidyverse)
library(data.table)

data_files <- dir(path = "./2020-07 IBS/data_225_frequent_output/")

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

#angles <- vector(mode = "list", length = length(data_files)) # use this if starting from scratch
#previous_files <- c() # use this if starting from scratch
angles <- readRDS("2020-07 IBS/angles_225_frequent_output_rds") # read in angles that have already been calculated since before
angles <- c(angles, vector(mode = "list", length = length(data_files)-length(angles))) # extend list so as to include new (not-yet-analysed) data files
previous_files <- angles %>% rbindlist(use.names = TRUE) %>% pull(file)
for(i in 1:length(data_files)){
   file <- data_files[i]
   
   if (file %notin% previous_files) {
      df <- fread(paste0("./2020-07 IBS/data_225_frequent_output/", file), col.names = c("t", "x", "y")) %>%   
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
            filter(diff >= 0) %>% ungroup() %>% # problem here: SOME REPLICATES HAVENT REACHED THIS!!!!! PUT SOME TEST HERE
            filter(t == min(t)) %>% pull(t) %>% unique()
         # test number of branches at end time, based on a sample of 1500 individuals
         db <- suppressMessages(df %>% 
                                   filter(t == branching_end_time) %>% 
                                   select(x, y) %>% 
                                   sample_n(1500) %>% 
                                   # We use DBSCAN: points within epsilon = ess/10 belong to same cluster; minimum points in a cluster is 100
                                   fpc::dbscan(ess/10, MinPts = 100, method = "raw"))
         
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

saveRDS(angles, file = "2020-07 IBS/angles_225_frequent_output_rds")


angles <- readRDS("2020-07 IBS/angles_225_frequent_output_rds")
output <- rbindlist(angles, use.names = TRUE) %>% as_tibble()

output %>% filter(nr_clusters == 2) %>% filter(onset) %>% nrow()

angle_onset_branching <- output %>% 
   filter(end_branching | onset) %>% 
   filter(nr_clusters == 2)

angle_before_onset <- filter(output, nr_clusters == 2)[which(filter(output, nr_clusters == 2)$onset)-1, ]$angle
angle_onset <- angle_onset_branching$angle[angle_onset_branching$onset]
angle_branching <- angle_onset_branching$angle[angle_onset_branching$end_branching]
plot(angle_before_onset, angle_branching)


exclude_from_distribution_data <- output[output$onset,] %>% 
   filter(t < 6000) %>% 
   pull(file)


distribution_data <- output %>% 
   filter(file %notin% exclude_from_distribution_data) %>% 
   filter(nr_clusters == 2) %>% 
   filter(t %in% c(0, 2000, 4000, 6000))

distribution_data <- rbind(distribution_data,
                           (distribution_data %>% mutate(angle = angle - pi)),
                            (distribution_data %>% mutate(angle = angle + pi))
                           )

distribution_data %>% 
   ggplot(aes(x = angle, y = t, group = t)) +
   ggridges::geom_density_ridges(alpha = 0.3, scale = 1.5) +
   geom_vline(xintercept = 3*pi/4) +
   coord_cartesian(xlim = c(0,pi)) + xlab("G matrix orientation (radians)") +
   ylab("Generations") +
   scale_y_continuous(expand = expansion(add = c(0, 500)),
                      breaks = c(0, 2000, 4000, 6000)) +
   scale_x_continuous(expand = c(0,0),
                      breaks = c(0, pi/4, pi/2, 3*pi/4, pi), 
                      labels = c(0, expression(frac(pi,4)), expression(frac(pi,2)), 
                                 expression(frac(3*pi,4)), expression(pi))) +
   cowplot::theme_cowplot() +
   theme(axis.text.y = element_text(vjust = 0))
