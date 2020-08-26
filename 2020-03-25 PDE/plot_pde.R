library(tidyverse)

ess <- 0.0813

c_x = 0.9; d_x = 1.65; b_x = 1.05;
b_y = b_x; c_y = c_x; d_y = d_x;	
ax = 1; ay = 1;
s = 100;

Bx <- function(x.mx, family) { 
  if (family == "power") { 
    benefits <- t((x.mx+t(x.mx))^alpha_x) 
    return(benefits)
  } else if (family == "quadratic") {
    benefits <- t((x.mx+t(x.mx))-b_x*(x.mx+t(x.mx))^2) 
    maxbx <- 1/(2*b_x) # coordinates of inflexion point
    benefits[x.mx+t(x.mx) >= maxbx] <- 1/(4*b_x)
    return(benefits)
  } else if (family == "mm") {
    benefits <- t(b1_x*(x.mx+t(x.mx))/(b2_x+(x.mx+t(x.mx))))
    return(benefits)
  }
}
By <- function(y.mx, family) { 
  if (family == "power") {
    benefits <- t((y.mx+t(y.mx))^alpha_y) 
    return(benefits)
  } else if (family == "quadratic") {
    benefits <- t((y.mx+t(y.mx))-b_y*(y.mx+t(y.mx))^2) 
    maxby <- 1/(2*b_y) # coordinates of inflexion point
    benefits[y.mx+t(y.mx) >= maxby] <- 1/(4*b_y)
    return(benefits) 
  } else if (family == "mm") {
    benefits <- t(b1_y*(y.mx+t(y.mx))/(b2_y+(y.mx+t(y.mx))))
    return(benefits)
  }
}
Cx <- function(x.ax, family) { 
  if (family == "power") {
    costs <- gamma_x*x.ax^beta_x     
    return(costs)
  } else if (family == "quadratic") {
    costs <- c_x * (x.ax - d_x * x.ax^2)
    maxcx <- 1/(2*d_x) # coordinates of inflexion point
    costs[x.ax >= maxcx] <- c_x/(4*d_x)
    return(costs)
  } else if (family == "mm") {
    costs <- c1_x*x.ax / (c2_x+x.ax)
    return(costs)
  }
}
Cy <- function(y.ax, family) { 
  if (family == "power") {
    costs <- gamma_y*y.ax^beta_y 
    return(costs)
  } else if (family == "quadratic") {
    costs <- c_y * (y.ax - d_y * y.ax^2) 
    maxcy <- 1/(2*d_y) # coordinates of inflexion point
    costs[y.ax >= maxcy] <- c_y/(4*d_y)
    return(costs)
  } else if (family == "mm") {
    costs <- c1_y*y.ax / (c2_y+y.ax)
    return(costs)
  }
}

directories <- dir("./2020-03-25 PDE/", pattern = paste0("zeta=2e-07_s=", s, "_mu=0.001"))

times <- c(0, 250, 350, 550)

angles_list <- vector(mode = "list", length = length(directories))

for(angle in 1:length(angles_list)){
  dir <- directories[angle]
  
  plot_list <- vector(mode = "list", length = length(times))
  
  for(i in 1:length(times)){
    timepoint <- readRDS(paste0("./2020-03-25 PDE/", dir, "/", times[i], "_full_population"))
    f <- timepoint[[1]]
    x.ax <- timepoint[[2]]
    y.ax <- timepoint[[3]]
    
    # Calculate fitness landscape
    grid <- length(x.ax)
    x.mx <- matrix(rep(x.ax,times=grid),grid,grid) # phenotype space matrix (x coordinate)
    y.mx <- matrix(rep(y.ax,each =grid),grid,grid) # (y coordinate)
    bx <- Bx(x.mx, "quadratic")
    cx <- Cx(x.ax, "quadratic")
    by <- By(y.mx, "quadratic")
    cy <- Cy(y.ax, "quadratic")
    xf <- apply(f,1,sum) # marginal distribution of population along x axis
    pix <- matrix(rep(apply(bx*xf, 2, sum) - cx, times = grid), grid, grid) # benefit - cost. bx*xf is frequency of each type of partner times their benefits.
    yf <- apply(f,2,sum)
    piy <- matrix(rep(apply(by*yf, 2, sum) - cy, each  = grid), grid, grid)
    pi <- ax*pix + ay*piy # payoff
    w <- s*pi; w[w<0] <- 0; # modified fdot
    mw <- sum(w*f) # mean fitness
    dw <- w-mw
    maxdw = max(dw)
    mindw <- min(dw)
    dw[dw<(-maxdw)] <- maxdw
    
    rownames(dw) <- x.ax
    colnames(dw) <- y.ax
    dw_long <- as.data.frame.table(dw, responseName = "fitness")
    colnames(dw_long) <- c("x", "y", "fitness")
    dw_long <- tibble(dw_long) %>% 
      mutate(x = as.numeric(as.character(x)), y = as.numeric(as.character(y)))
    
    rownames(f) <- x.ax
    colnames(f) <- y.ax
    f_long <- as.data.frame.table(f, responseName = "frequency")
    colnames(f_long) <- c("x", "y", "frequency")
    f_long <- tibble(f_long) %>% 
      mutate(x = as.numeric(as.character(x)), y = as.numeric(as.character(y)))
    
    
    dw_long <- dw_long %>% mutate(fitness2 = if_else(fitness > 0, 1, 0))
    
    # Plot
    ggplot() +
      geom_raster(data = dw_long, aes(x = x, y = y, fill = fitness)) +
      scale_fill_gradientn(
        colors = c("tomato","white","steelblue"),
        values = scales::rescale(c(mindw,0,maxdw)),
        limits = c(mindw,maxdw)
      ) +
      geom_contour(data = f_long, aes(x = x, y = y, z = frequency), color = "black", bins = 6) +
      cowplot::theme_cowplot() +
      theme(legend.position = "none",
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 9),
            aspect.ratio = 1) +
      scale_x_continuous(expand = expansion(0,0),
                         labels = scales::label_number(accuracy = 0.01)) +
      scale_y_continuous(expand = expansion(0,0),
                         labels = scales::label_number(accuracy = 0.01)) +
      geom_hline(yintercept = ess) +
      geom_vline(xintercept = ess) +
      xlab(NULL) + ylab(NULL) -> plot
    
    if (i == 1) {
      plot <- plot + theme(axis.title.y = element_text(size = 11)) + labs(y = expression(z[2]))
      if(angle == 1) plot <- plot + labs(tag = "A")
      if(angle == 2) plot <- plot + labs(tag = "B")
      if(angle == 3) plot <- plot + labs(tag = "C")
      if(angle == 4) plot <- plot + labs(tag = "D")
      }
    if(i > 1) plot <- plot + theme(axis.text.y = element_blank())
    if(angle == 4) plot <- plot + theme(axis.title.x = element_text(size = 11)) + labs(x = expression(z[1]))
    
    
    plot_list[[i]] <- plot
  }
  angles_list[[angle]] <- plot_list
}


plot_pde <- egg::ggarrange(angles_list[[1]][[1]], angles_list[[1]][[2]], 
               angles_list[[1]][[3]], angles_list[[1]][[4]], 
               angles_list[[2]][[1]], angles_list[[2]][[2]], 
               angles_list[[2]][[3]], angles_list[[2]][[4]],
               angles_list[[3]][[1]], angles_list[[3]][[2]], 
               angles_list[[3]][[3]], angles_list[[3]][[4]],
               angles_list[[4]][[1]], angles_list[[4]][[2]], 
               angles_list[[4]][[3]], angles_list[[4]][[4]],
               nrow = 4, ncol = 4,
               labels = c("A", "", "", "",
                          "B", "", "", "",
                          "C", "", "", "",
                          "D", "", "", "")
)


ggsave(filename = "./2020-03-25 PDE/plot_pde_2.pdf",
       plot = plot_pde,
       width = 8.5, height = 8.5)
