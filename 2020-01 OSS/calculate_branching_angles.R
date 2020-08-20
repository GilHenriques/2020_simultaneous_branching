library(tidyverse)
library(data.table)

ess <- 0.0813
nr_poly <- 10 # nr of invasions with polymorphic population until branching is "complete"

output_from_OSS <- dir(path = "./2020-01 OSS/data/", pattern = "output_list")
angles_of_approach <- parse_number(output_from_OSS)

branching_angles_list <- vector(mode = "list", length = length(angles_of_approach))

for(i in seq_along(output_from_OSS)){
  data <- paste0("./2020-01 OSS/data/", output_from_OSS[i]) %>% 
    read_rds() %>% 
    map(. %>% drop_na())
  
  branching_angles_df <- tibble(angle_of_approach = angles_of_approach[i], rep = seq_along(data), angle_of_branching = 999)
  
  for(rep in seq_along(data)){
    time_strains <- suppressMessages(data[[rep]] %>% 
                                       drop_na() %>% 
                                       group_by(t) %>% 
                                       summarize(nr_strains = n()) %>% 
                                       mutate(poly = nr_strains >= 2) 
    )
    
    counter <- 0; time_strains$sequential_poly <- 0;
    for(j in 1:nrow(time_strains)){
      if(time_strains[j, ]$poly) { counter <- counter + 1 } else { counter <- 0 }
      time_strains[j, ]$sequential_poly <- counter
    }
    
    time_branching <- time_strains %>% 
      filter(sequential_poly == nr_poly) %>% 
      slice(1) 
    
    if(nrow(time_branching) == 0 || time_branching$nr_strains != 2) {
      branching_angles_df[rep, ]$angle_of_branching <- NA
    } else {
      time_branching <- time_branching %>% pull(t)
      
      points <- data[[rep]] %>% filter(t == time_branching) %>% select(x,y) 
      
      angle <- atan2((points[1, ]$y - points[2, ]$y), (points[1, ]$x - points[2, ]$x))  
      if (angle < 0) angle <- angle + pi
      
      branching_angles_df[rep, ]$angle_of_branching <- angle
    }
    if(rep %% 50 == 0) print( paste0(round(rep/length(data)*100, 2), "% of treatment ", i, "/", length(output_from_OSS), " complete!"))
  }
  
  branching_angles_list[[i]] <- branching_angles_df
}

all_branching_angles <- bind_rows(branching_angles_list)
saveRDS(all_branching_angles, file = "./2020-01 OSS/branching_angles_NEW_rds")
