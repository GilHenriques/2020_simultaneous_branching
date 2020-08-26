# Asymmetric games

library(tidyverse)
library(raster)

directories <- dir(path = "2020-07 PDE asymmetry/data_for_comparing_with_ibs/", pattern = "2020-07-29")
threshold <- 1e-6
multimodal <- function(f){ # returns true if the marginal distribution is multimodal along either x or y
  f[f < threshold] <- 0
  marginal_y <- colSums(f); marginal_x <- rowSums(f)
  multimodal_y <- length(rle(marginal_y > (max(marginal_y))/10)$values) > 3
  multimodal_x <- length(rle(marginal_x > (max(marginal_x))/10)$values) > 3
  return(any(multimodal_x, multimodal_y))
}

angles_df <- tibble(dir = directories, angle = 9999)
for(directory in directories){
  timepoints <- dir(path = paste0("2020-07 PDE asymmetry/data_for_comparing_with_ibs/", directory), pattern = "full_population")
  timepoints <- gtools::mixedsort(timepoints)
  maxXY <- NA
  
  for (timepoint in timepoints){
    df <- readRDS(paste("2020-07 PDE asymmetry/data_for_comparing_with_ibs", directory, timepoint, sep = "/"))[[1]]
    if (multimodal(df)) { 
      r <- raster(t(df)) %>% flip(direction = "y")
      extent(r) <- extent(c(0, nrow(df), 0, nrow(df)) + 0.5)
      
      ## Find the maximum value within the 9-cell neighborhood of each cell
      f <- function(X) max(X, na.rm=TRUE)
      ww <- matrix(1, nrow=3, ncol=3) ## Weight matrix for cells in moving window
      localmax <- focal(r, fun=f, w=ww, pad=TRUE, padValue=NA)
      
      ## Does each cell have the maximum value in its neighborhood?
      r2 <- r==localmax
      
      ## Get x-y coordinates of those cells that are local maxima
      maxXY <- xyFromCell(r2, Which(r2==1 & r != 0, cells=TRUE))
      
      print(timepoint)
      break
    }
  }
  maxXY <- maxXY %>% as_tibble()
  x1 <- maxXY[1,]$x; y1 <- maxXY[1,]$y; x2 <- maxXY[2,]$x; y2 <- maxXY[2,]$y
  varx <- var(c(x1, x2)); vary <- var(c(y1,y2)); cov <- cov(c(x1,x2),c(y1,y2));
  evec_x = - (- varx + vary - sqrt(4*cov^2 + varx^2 - 2*varx*vary + vary^2 ))/(2*cov)
  if(!is.finite(evec_x)) evec_x <- 0
  evec_y = 1
  angles_df[angles_df$dir == directory,]$angle <- acos(pracma::dot(c(evec_x, evec_y), c(1,0))/(pracma::Norm(c(evec_x, evec_y))*pracma::Norm(c(1,0))))
}

angles_df <- angles_df %>% 
  separate(dir, into = c("dir", "alpha"), sep = "alpha=") 
angles_df[angles_df$alpha == 0, ]$angle <- pi # fix numerical mistake

saveRDS(angles_df, "2020-07 PDE asymmetry/angle_vs_y_for_comparing_with_ibs_rds")

