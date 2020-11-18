library(tidyverse)
library(data.table)

pop_size <- 10000
bx <- 1.05; cx <- 0.9; dx <- 1.65; ess <-  - (cx - 1)/(2 * (2*bx - cx*dx))

angles <- list.files("2020-10_IBS_endstate/data")

for (j in seq_along(angles)) {
   data_files <- dir(path = paste0("./2020-10_IBS_endstate/data/", angles[j]))

   output <- vector(mode = "list", length = length(data_files))
   
   for(i in 1:length(data_files)){
      file <- data_files[i]
      
      df <- fread(paste0("./2020-10_IBS_endstate/data/", angles[j], "/", file))
      
      if(ncol(df) == 3) {
         colnames(df) <- c("t", "x", "y")
         
         df <- df %>% filter(t == max(t)) %>% 
            mutate(t = t/pop_size)  # convert time to generations
         
         # test number of branches at end time, based on a sample of 1500 individuals
         db <- suppressMessages(df %>% 
                                   select(x, y) %>% 
                                   sample_n(1500) %>% 
                                   # We use DBSCAN: points within epsilon = ess/10 belong to same cluster; minimum points in a cluster is 100
                                   fpc::dbscan(ess/10, MinPts = 100, method = "raw"))
         
         nr_clusters <- db$cluster %>% unique() %>% length()
         
         output[[i]]$file <- file
         output[[i]]$nr_clusters <- nr_clusters
         
         if(nr_clusters < 2) {
            output[[i]]$class <- "Other"
         } else {
            
            null <- FALSE; x <- FALSE; y <- FALSE; xy <- FALSE; 
            if (filter(df, x > ess & y > ess) %>% nrow() > 0) xy <- TRUE
            if (filter(df, x < ess & y < ess) %>% nrow() > 0) null <- TRUE
            if (filter(df, x > ess & y < ess) %>% nrow() > 0) x <- TRUE
            if (filter(df, x < ess & y > ess) %>% nrow() > 0) y <- TRUE
            
            output[[i]]$class <- if(null & xy) { "Cooperator & defector (CD)"
            } else if (x & y) { "Division of labor (DL)" 
            } else if (x & y & null) { "DL + Defector"
            } else if (x & y & xy) { "DL + Cooperator"
            } else if (null & x & xy) {"CD + Specialist 1"
            } else if (null & y & xy) { "CD + Specialist 2"
            } else if (null & x & y & xy) { "Four types" }
            
         }
         
         output[[i]]$angle_approach <- angles[j]
      }
      
      print(paste0("folder ", angles[j], " file ", i))
   }
   saveRDS(output, file = paste0("2020-10_IBS_endstate/endstates_", angles[j], "_rds"))

}



