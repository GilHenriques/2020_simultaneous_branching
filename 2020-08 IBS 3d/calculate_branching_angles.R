library(tidyverse)
library(data.table)

data_files <- dir(path = "./2020-08 IBS 3d/data")

b <- 1.05; c <- 0.9; d <- 1.65; ess <-  - (c - 1)/(2 * (2*b - c*d))

x0 <- 0.060087; y0 <- 0.060087; z0 <- 0.060087

# Vector N defines the direction of approach
N <- c(x0 - ess, y0 - ess, z0 - ess)


# Plane Q includes the arbitrary point A = (1,2,3), the initial position (x0, y0, z0), and the ESS.
# Plane Q is defined by its normal vector M.
M <- pracma::cross(N, c(1-ess, 2-ess, 3-ess))

pop_size <- 10000
`%notin%` <- Negate(`%in%`)

get_fitness <- function(x, y, xbar, ybar){
  bx <- (x + xbar) - b * (x + xbar)^2
  by <- (y + ybar) - b * (y + ybar)^2
  cx <- c * (x - d * x^2)
  cy <- c * (y - d * y^2)
  return((bx-cx) + (by-cy))
}


angles <- vector(mode = "list", length = length(data_files)) # use this if starting from scratch
previous_files <- c()
for(i in 1:length(data_files)){
  file <- data_files[i]
  
  if (file %notin% previous_files) {
    df <- fread(paste0("./2020-08 IBS 3d/data/", file), col.names = c("t", "x", "y", "z")) %>%   
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
      
      df_sample <- suppressMessages(df %>% filter(t == branching_end_time) %>% 
        select(x, y, z) %>% 
        sample_n(1500))
      
      db <- suppressMessages(df_sample %>% 
                               # We use DBSCAN: points within epsilon = ess/10 belong to same cluster; minimum points in a cluster is 100
                               fpc::dbscan(ess/10, MinPts = 100, method = "raw"))
    
      nr_clusters <- db$cluster %>% unique() %>% length()
      
      angles[[i]]$file <- file
      angles[[i]]$nr_clusters <- nr_clusters
      
      if(nr_clusters != 2) {
        angles[[i]]$t <- NA
        angles[[i]]$alpha <- NA
        angles[[i]]$beta <- NA
      } else {
        df_sample$cluster <- db$cluster
       
        pop <- suppressMessages(df_sample %>% group_by(cluster) %>% summarize(mean_x = mean(x), mean_y = mean(y), mean_z = mean(z)) )
        
        x1 <- pop[1,]$mean_x; y1 <- pop[1,]$mean_y; z1 <- pop[1,]$mean_z
        x2 <- pop[2,]$mean_x; y2 <- pop[2,]$mean_y; z2 <- pop[2,]$mean_z
        
        # vector L: defined by the two branches
        L <- c(x1 - x2, y1 - y2, z1 - z2)
        # Vector V: defined by the intersection of P and Q, i.e., V = N x M (where x denotes cross product)
        V <- pracma::cross(N, M)
        
        # The vector LP is the projection of L on P. 
        # It's defined by vector N and vector L, i.e. LP = L - (L.N/(||N||^2)) * N (where . is the dot product)
        LP <- L - (pracma::dot(L, N)/(pracma::Norm(N)^2))*N
        
        # The vector LQ is the projection of L on Q. 
        # It's defined by vector Q and vector L, i.e. LQ = L - (L.M/(||M||^2)) * M (where . is the dot product)
        LQ <- L - (pracma::dot(L, M)/(pracma::Norm(M)^2))*M
        
        # Alpha is the angle between LQ and V. cos alpha = LQ.V/(||LQ||*||V||)
        # Beta is the angle between LP and V. cos beta = LP.V/(||LP||*||V||)
        angles[[i]]$t <- branching_end_time
        angles[[i]]$alpha <- (pracma::dot(LQ, V)/(pracma::Norm(LQ)*pracma::Norm(V))) %>% acos
        angles[[i]]$beta <- (pracma::dot(LP, V)/(pracma::Norm(LP)*pracma::Norm(V))) %>% acos
        }
    } else {
      # nothing happens 
      angles[[i]]$file <- file
      angles[[i]]$nr_clusters <- NA
      angles[[i]]$t <- NA
      angles[[i]]$alpha <- NA
      angles[[i]]$beta <- NA
    }
    print(i)
  }
}

saveRDS(angles, file = "2020-08 IBS 3d/angles_rds")
