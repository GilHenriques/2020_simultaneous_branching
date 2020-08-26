# Packages ----------------------------------------------------------------
if(!require("tidyverse")) install.packages("tidyverse"); require("tidyverse")

# Bivariate normal --------------------------------------------------------
dnorm2 <- function(x, y, ux = 0, uy = 0, sx = 1, sy = 1, rho = 0) {
  1/(2*acos(-1)*sx*sy*sqrt(1-rho*rho))*exp(-1/(2*(1-rho*rho))*((x-ux)^2/sx/sx - 2*rho*(x-ux)*(y-uy)/sx/sy+(y-uy)^2/sy/sy))
}

# Replicator-mutator equation ---------------------------------------------
fdot <- function(f, w, mutation_rate, grid, y.ax, x.ax, mw, zeta, bs) {
  # part 1: Sum_j (f_j * w_j * Mutation from i to j)
  change <- f*w
  change <- (1-(mutation_rate/bs^2))*change +
	  mutation_rate/4/bs^2*(  # used 4 (number of neighbours) for von Neumann neighborhood
	    rbind(change[-1,],matrix(0,1,grid))+
      rbind(matrix(0,1,grid),change[-grid,])+
      cbind(change[,-1],matrix(0,grid,1))+
      cbind(matrix(0,grid,1),change[,-grid]) +
	    (y.ax[1] < 10e-10) * cbind(change[,1], matrix(0,grid,grid-1)) +
      (x.ax[1] < 10e-10) * rbind(change[1,], matrix(0,grid-1,grid))
  )

  # part 2: (wbar * f_i) - zeta (zeta is to make sure population doesn't take up all the pht space)
  change <- change - mw*f - zeta/bs^2
  return(change)
}


# Summarize ---------------------------------------------------------------
summarize <- function(x.ax, y.ax, grid, f, density_summary){

  max.row.col <- as.data.frame(matrix(NA,ncol=2,nrow=0))

  # find coordinates of maxima
  for (row in 1:grid){
    for (col in 1:grid) {
      nr.neighbours <- 0; counter <- 0;
      if (row < grid) {
        nr.neighbours <- nr.neighbours+1
        if (f[row,col] > f[row+1,col]) counter <- counter+1
      }
      if (row > 1) {
        nr.neighbours <- nr.neighbours+1
        if (f[row,col] > f[row-1,col]) counter <- counter+1
      }
      if (col < grid) {
        nr.neighbours <- nr.neighbours+1
        if (f[row,col] > f[row,col+1]) counter <- counter+1
      }
      if (col > 1) {
        nr.neighbours <- nr.neighbours+1
        if (f[row,col] > f[row,col-1]) counter <- counter+1
      }
      if (nr.neighbours == counter) { # local maximum
        max.row.col <- rbind(max.row.col, c(row, col))
      }
    }
  }

  # create lookup matrix
  lookup <- matrix(NA, nrow = grid, ncol = grid)
  for (i in 1:nrow(max.row.col)){
    lookup[max.row.col[i,1],max.row.col[i,2]] <- i
  }

  summary <- as.data.frame(which(lookup < 100, arr.ind = T))
  summary$id <- NA
  for (i in 1:nrow(summary)){
    summary[i,3] <- lookup[summary[i,1],summary[i,2]]
  }
  summary <- data.frame(id = summary$id, y = y.ax[summary[,2]], x = x.ax[summary[,1]])
  summary$density <- NA

  if (density_summary == TRUE){
    # calculate densitites
    thresh <- 1e-10

    while(sum(lookup < 100, na.rm=T) > 0) {
      to.check <- which(lookup < 100, arr.ind = T)
      for (check.point in 1:nrow(to.check)){
        count <- 0
        if(to.check[check.point,1]+1 <= grid) {
          if (
            f[to.check[check.point,1],to.check[check.point,2]] >= f[to.check[check.point,1]+1,to.check[check.point,2]] &
            f[to.check[check.point,1]+1,to.check[check.point,2]] > thresh
          ) {
            lookup[to.check[check.point,1]+1,to.check[check.point,2]] <- lookup[to.check[check.point,1],to.check[check.point,2]]
          }
        }
        count <- count+1

        if(to.check[check.point,1]-1 >= 1) {
          if (
            f[to.check[check.point,1],to.check[check.point,2]] >= f[to.check[check.point,1]-1,to.check[check.point,2]] &
            f[to.check[check.point,1]-1,to.check[check.point,2]] > thresh
          ) {
            lookup[to.check[check.point,1]-1,to.check[check.point,2]] <- lookup[to.check[check.point,1],to.check[check.point,2]]
          }
        }
        count <- count+1

        if(to.check[check.point,2]+1 <= grid) {
          if(
            f[to.check[check.point,1],to.check[check.point,2]] >= f[to.check[check.point,1],to.check[check.point,2]+1] &
            f[to.check[check.point,1],to.check[check.point,2]+1] > thresh
          ) {
            lookup[to.check[check.point,1],to.check[check.point,2]+1] <- lookup[to.check[check.point,1],to.check[check.point,2]]
          }
        }
        count <- count+1

        if(to.check[check.point,2]-1 >= 1) {
          if(
            f[to.check[check.point,1],to.check[check.point,2]] >= f[to.check[check.point,1],to.check[check.point,2]-1] &
            f[to.check[check.point,1],to.check[check.point,2]-1] > thresh
          ) {
            lookup[to.check[check.point,1],to.check[check.point,2]-1] <- lookup[to.check[check.point,1],to.check[check.point,2]]
          }
        }
        count <- count+1

        if (count == 4) { lookup[to.check[check.point,1],to.check[check.point,2]] <- lookup[to.check[check.point,1],to.check[check.point,2]]*100 }
      }
    }

    for (i in 1:nrow(summary)){
      code <- summary[i,1]*100
      summary[i,4] <- sum(f[lookup == code],na.rm = T)
    }
  } # end density summary

  return(summary)
}


# Costs and benefits ------------------------------------------------------
Bx <- function(x.mx, family) {
  if (family == "power") { # didn't use
    benefits <- t((x.mx+t(x.mx))^alpha_x)
    return(benefits)
  } else if (family == "quadratic") {
    benefits <- t((x.mx+t(x.mx))-b_x*(x.mx+t(x.mx))^2)
    maxbx <- 1/(2*b_x) # coordinates of inflexion point
    benefits[x.mx+t(x.mx) >= maxbx] <- 1/(4*b_x)
    return(benefits)
  } else if (family == "mm") { # didn't use
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


# ESS ---------------------------------------------------------------------
ess_x <- function(family){
  if (family == "power") { ((alpha_x *2^(alpha_x-1))/(gamma_x*beta_x))^(1/(beta_x - alpha_x))
  } else if (family == "quadratic") { - (c_x - 1)/(2*(2*b_x-c_x*d_x))
  } else if (family == "mm") {
    if (2*c2_x > b2_x) {
      (-b2_x*(b1_x-2*c1_x)*c2_x - sqrt(b1_x*b2_x*c1_x*c2_x*(b2_x-2*c2_x)^2))/(b1_x*b2_x-4*c1_x*c2_x)
    } else if (2*c2_x < b2_x){
      (-b2_x*(b1_x-2*c1_x)*c2_x + sqrt(b1_x*b2_x*c1_x*c2_x*(b2_x-2*c2_x)^2))/(b1_x*b2_x-4*c1_x*c2_x)
    }
  }
}
ess_y <- function(family){
  if (family == "power") { ((alpha_y *2^(alpha_y-1))/(gamma_y*beta_y))^(1/(beta_y - alpha_y))
  } else if (family == "quadratic") { - (c_y - 1)/(2*(2*b_y-c_y*d_y))
  } else if (family == "mm") {
    if (2*c2_y > b2_y) {
      (-b2_y*(b1_y-2*c1_y)*c2_y - sqrt(b1_y*b2_y*c1_y*c2_y*(b2_y-2*c2_y)^2))/(b1_y*b2_y-4*c1_y*c2_y)
    } else if (2*c2_y < b2_y){
      (-b2_y*(b1_y-2*c1_y)*c2_y + sqrt(b1_y*b2_y*c1_y*c2_y*(b2_y-2*c2_y)^2))/(b1_y*b2_y-4*c1_y*c2_y)
    }
  }
}


# PDE function ------------------------------------------------------------
play.evol = function(
  Angle = 225, # angle of origin. 225 = bottom left
  mutation_rate0 = 1e-8, s = 80, zeta0 = 2e-14, # s = strength of selection; zeta controls how much the population can "spread" around fitness maxima
  abserr = 1e-5,   # permissible absolute error, i.e., |x-y| (for runge-kutta calculation)
  relerr = 1e-4,   # permissible relative error, i.e., |x-y|/|x| (for runge-kutta calculation)
  err_order = 4,   # constant: runge_kutta dopri5's error order is 4
  step_order = 5,  # constant: runge_kutta dopri5's step order is 5
  dt0 = 0.001,     # initial value of dt (for runge kutta calculation)
  threshold = 1e-6, # threshold for detecting stable state - Stable state detection doesn't work very well
  unstable_ratio = 1e+1, unstable_threshold = 2e-3, # used for detecting stable state when time step is fluctuating
  itstep = 100, itermax = 300*itstep,  # print output every itstep iterations (one iteration = dt, but dt changes dynamically)
  rho = 0, asym = 0.5, # multiply fitness of each game
  figures = TRUE, # save figures every few iterations if TRUE, else only final time-step
  density_summary = FALSE, # save x and y density if TRUE, else doesn't save densities
  family = "quadratic", # can be power or quadratic or mm
  epistasis = "additive", # can be additive or multiplicative
  name = "test") {

  # Set up output directory and output matrix -------------------------------
  dir.name <- name; dir.create(file.path(getwd(), dir.name))

  # output matrix will be 10 times as long as the number of output events to avoid running out of space
  output <- matrix(NA, ncol = 8, nrow = 10*ceiling(itermax/itstep)+1)
  colnames(output) <- c("id", "Y", "X", "density", "t", "dt", "err1", "err2")

  # parameters
  xstar <- ess_x(family); ystar <- ess_y(family)
  zeta <- zeta0*(xstar/100)^2
  mutation_rate <- mutation_rate0*(xstar/100)^2

  # Set up grid -------------------------------------------------------------
  grid <- 0; bs <- xstar/100; # phenotype distribution bin size
  while(grid < 125 | grid > 250){
    dist <- (min(xstar,ystar))/2 # distance between initial population and EBP
    #dist <- 0.03 # distance between initial population and EBP, used this value for supplemental PDE with correlation
    x0 <- xstar + dist * cos(Angle/180 * acos(-1)) # initial population
    y0 <- ystar + dist * sin(Angle/180 * acos(-1))

    # window will be of width 3*dist, if possible centred around x0. Rounding is so that 0 is a possible value for the phenotype
    min.x <- max(0, ceiling((x0 - (dist+dist/2))/bs)*bs)
    max.x <- min.x + ceiling(3*dist/bs)*bs
    min.y <- max(0, ceiling((y0 - (dist+dist/2))/bs)*bs)
    max.y <- min.y + ceiling(3*dist/bs)*bs

    x.ax <- seq(min.x, max.x, by = bs) # phenotype axes
    y.ax <- seq(min.y, max.y, by = bs)

    grid <- length(x.ax) # should be equal to length of y.ax
    if(grid%%2 != 0){ # must be possible to divide grid in half
      max.x <- max.x + bs
      max.y <- max.y + bs
      x.ax <- seq(min.x, max.x, by = bs) # phenotype axes
      y.ax <- seq(min.y, max.y, by = bs)
      grid <- length(x.ax)
    }

    if (grid < 125){ bs <- bs/2
    } else if (grid > 250) { bs <- 2*bs }
  }

  x.mx <- matrix(rep(x.ax,times=grid),grid,grid) # phenotype space matrix (x coordinate)
  y.mx <- matrix(rep(y.ax,each =grid),grid,grid) # (y coordinate)


  # Set up initial population -----------------------------------------------
  pi <- matrix(0, grid, grid) # matrix that will hold payoffs of phenotypes
  w <- matrix(0, grid, grid) # matrix that will hold fitnesses of phenotypes
  f <- matrix(0, grid, grid) # will hold distribution of individuals
  f <- dnorm2(x.mx, y.mx, x0, y0, 0.01*xstar, 0.01*ystar, rho) # start w normal distribution of individuals
  #f <- dnorm2(x.mx, y.mx, x0, y0, 0.05*xstar, 0.05*ystar, 0.9) # this was used for starting with positive rho
  f <- f/sum(f) # integral is one

  # these will be used to check whether we should interrupt simulation
  f_old <- matrix(0, grid, grid) # Will hold f(iter-itstep)
  f_older <- matrix(1, grid, grid) # Will hold f(iter-2*itstep)


  # Costs and benefits ------------------------------------------------------
  bx <- Bx(x.mx, family)
  cx <- Cx(x.ax, family)
  by <- By(y.mx, family)
  cy <- Cy(y.ax, family)

  # Begin iteration ---------------------------------------------------------
  itercounter <- itstep; # counts from zero to itstep, and print output whenever itercounter = itstep
  iter <- 0; # counts from zero to itermax
  t <- 0;  dt <- dt0;

  while(iter <= itermax) {

    # Expand matrix if needed -------------------------------------------------
    if (sum(f[c(1,grid),])+sum(f[,c(1,grid)]) > 0){ # if the distribution reaches the edges of phenotype space

      if ( sum(f[grid,]) > 0 | sum(f[,grid]) > 0 ){ # population reached right or top edge
        f  <- f[c(T,F),c(T,F)] + f[c(F,T),c(F,T)] + f[c(T,F), c(F,T)] + f[c(F,T),c(T,F)] # f is now half as big

        f <- rbind(f, matrix(0, nrow = grid/2, ncol = ncol(f)))
        f <- cbind(f, matrix(0, nrow = nrow(f), ncol = grid/2))

        bs <- bs*2
        x.ax <- seq(from = x.ax[1], length = grid, by = bs)
        y.ax <- seq(from = y.ax[1], length = grid, by = bs)

        x.mx <- matrix(rep(x.ax,times=grid),grid,grid) # phenotype space matrix (x coordinate)
        y.mx <- matrix(rep(y.ax,each =grid),grid,grid) # (y coordinate)

        # Costs and benefits
        bx <- Bx(x.mx, family)
        cx <- Cx(x.ax, family)
        by <- By(y.mx, family)
        cy <- Cy(y.ax, family)
      }


      if ((x.ax[1] > 0 & sum(f[1,]) > 0) | (y.ax[1] > 0 & sum(f[,1]) > 0)){ # population reached left edge or bottom (increase x axis but no more than zero)
        startx <- max(0, x.ax[1]-bs*grid)
        starty <- max(0, y.ax[1]-bs*grid)

        bs <- bs*2
        if (startx == 0 | starty == 0) {
          if (x.ax[1]%%bs != 0 | y.ax[1]%%bs != 0) {
            f = rbind(rep(0, length = ncol(f)), f, rep(0, length = ncol(f)))
            f = cbind(rep(0, length = nrow(f)), f, rep(0, length = nrow(f)))
            x.ax = c(x.ax-bs/2, x.ax, tail(x.ax,1)+bs/2)
            y.ax = c(y.ax-bs/2, y.ax, tail(y.ax,1)+bs/2)
          }
        }

        f <- f[c(T,F),c(T,F)] + f[c(F,T),c(F,T)] + f[c(T,F), c(F,T)] + f[c(F,T),c(T,F)] # f is now half as big

        if(x.ax[1]-bs <=0) { newlowx <- 0
        } else { newlowx <- length(seq(from = startx, to = x.ax[1]-bs, by = bs)) }
        if(y.ax[1]-bs <=0) { newlowy <- 0
        } else { newlowy <- length(seq(from = starty, to = y.ax[1]-bs, by = bs)) }
        newhighx <- grid - (ncol(f) + newlowx)
        newhighy <- grid - (nrow(f) + newlowy)

        f <- rbind(matrix(0, nrow = newlowx, ncol = ncol(f)), f, matrix(0, nrow = newhighx, ncol = ncol(f)))
        f <- cbind(matrix(0, nrow = nrow(f), ncol = newlowy), f, matrix(0, nrow = nrow(f), ncol = newhighy))

        x.ax <- seq(from = startx, by = bs, length = grid)
        y.ax <- seq(from = starty, by = bs, length = grid)

        x.mx <- matrix(rep(x.ax,times=grid),grid,grid) # phenotype space matrix (x coordinate)
        y.mx <- matrix(rep(y.ax,each =grid),grid,grid) # (y coordinate)

        # Costs and benefits
        bx <- Bx(x.mx, family)
        cx <- Cx(x.ax, family)
        by <- By(y.mx, family)
        cy <- Cy(y.ax, family)
      }
    }

    # Calculate payoff and fitness --------------------------------------------
    xf <- apply(f,1,sum) # marginal distribution of population along x axis
    pix <- matrix(rep(apply(bx*xf, 2, sum) - cx, times = grid), grid, grid) # benefit - cost. bx*xf is frequency of each type of partner times their benefits.

    yf <- apply(f,2,sum)
    piy <- matrix(rep(apply(by*yf, 2, sum) - cy, each  = grid), grid, grid)

    pi <- (1-asym)*pix + (asym)*piy # payoff
    if(epistasis == "multiplicative") pi <- (1-asym)*pix * (asym)*piy # payoff
    w <- s*pi; w[w<0] <- 0; # modified fdot
    mw <- sum(w*f) # mean fitness

    # Printing ----------------------------------------------------------------
    if (itercounter == itstep){
      print(t)

      # .  Reset counter and check if stop simulation ------------------------------
      itercounter  <- 0

      err1 <- max(abs(f-f_old))
      err2 <- max(abs(f-f_older))
      if (
        (iter/itstep+1 > 30) && # don't end simulation before 30 snapshots have been produced
        (err1 < threshold && err2 < threshold) || (err1/err2 > unstable_ratio && err1 < unstable_threshold)
      ) { iter <- itermax }

      f_older <- f_old # two steps ago becomes one step ago
      f_old <- f # one step ago becomes now

      # .  Calculate summary -------------------------------------------------------
      summr <- summarize(x.ax, y.ax, grid, f, density_summary)
      summr$t <- t
      summr$dt <- dt
      summr <- as.matrix(summr)
      # save err1, err2 to determine threshold
      non.na.rows <- rowSums(is.na(output)) != ncol(output)
      first.row <- which(non.na.rows == F)[1] # first empty row
      output[first.row:(first.row+nrow(summr)-1), ncol(output)-1] <- err1
      output[first.row:(first.row+nrow(summr)-1), ncol(output)] <- err2
      output[first.row:(first.row+nrow(summr)-1), 1:6] <- summr

      # .  Make figures ------------------------------------------------------------
      if (figures == TRUE){ #| iter == itermax) {
        dw = w-mw
        maxdw = max(w-mw)
        dw[dw<(-maxdw)]=-maxdw

        png(paste("./", dir.name, "/", iter, ".png", sep=""), height = 1600, width = 1600)
        par(cex = 2.0)
        image(x.ax,y.ax,dw,zlim=c(-maxdw, maxdw), col = cm.colors(12), main = t)
        contour(x.ax, y.ax, f, add = TRUE)
        draw_f = f; draw_f[draw_f > 0] <- 1
        contour(x.ax, y.ax, draw_f, col="green", add=TRUE)
        abline(v = xstar, lty = 2)
        abline(h = xstar, lty = 2)
        dev.off()

		saveRDS(list(f,x.ax,y.ax,t,dt), paste("./", dir.name,"/", iter, "_full_population",sep=""))
      }
    }

    # Calculate new frequency (Runge-Kutta-Fehlberg) --------------------------
    while(TRUE){
      # Part 1: calculate the five slopes of the Runge-Kutta-Fehlberg method
      k1 <- fdot(f, w, mutation_rate, grid, y.ax, x.ax, mw, zeta, bs)
      k2 <- fdot(f + dt/4         * k1, w, mutation_rate, grid, y.ax, x.ax, mw, zeta, bs)
      k3 <- fdot(f + 3*dt/32      * k1 + 9*dt/32      * k2, w, mutation_rate, grid, y.ax, x.ax, mw, zeta, bs)
      k4 <- fdot(f + 1932*dt/2197 * k1 - 7200*dt/2194 * k2 + 7296*dt/2197 * k3, w, mutation_rate, grid, y.ax, x.ax, mw, zeta, bs)
      k5 <- fdot(f + 439*dt/216   * k1 - 8*dt         * k2 + 3680*dt/513  * k3 - 845*dt/4104  * k4, w, mutation_rate, grid, y.ax, x.ax, mw, zeta, bs)
      k6 <- fdot(f - 8*dt/27      * k1 + 2*dt         * k2 - 3544*dt/2565 * k3 + 1859*dt/4104 * k4 - 11*dt/40 * k5, w, mutation_rate, grid, y.ax, x.ax, mw, zeta, bs)

      # Part 2: calculate the candidates for f
      f4order <- f + dt*(25*k1/216 + 1408*k3/2565  + 2197*k4/4104   - k5/5)
      f5order <- f + dt*(16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55)

      # Part 3: calculate the local truncation error
      err   <- max(abs(f4order - f5order)/(abserr+relerr*(abs(f)+abs(f4order-f))))

      # Part 4: update time and adaptively update dt
      if(err < 1){
        # Success! All error elements of are less than threshold
        t  <- t+dt

        if(err < 0.5){
          # Try increase time step
          err  <- max(5^(-step_order),err)
          dt <- dt*(9/10*err^(-1/step_order))
        }
        break
      }
      # Fail! Some error elements are larger than threshold
      # Decrease time step
      dt <- dt*max(9/10*err^(-1/(err_order-1)),1/5)

      # Retry Runge-Kutte calculation
    }


    # Final updates -----------------------------------------------------------
    f <- f4order
    f[f<0] <- 0
    f <- f/sum(f)

    itercounter <- itercounter+1; # itercounter goes from 0 to itstep, then resets
    iter <- iter+1; # iter goes from 0 to itermax

  } # end iteration


  # Save parameters and RDS -------------------------------------------------
  params <- data.frame(
    name = c("epistasis", "angle", "mutation_rate0", "s", "zeta0", "abserr", "relerr", "err_order", "step_order", "dt0", "itstep", "itermax", "threshold", "unstable_ratio", "unstable_threshold", "rho", "asym", "figures", "family", "density_summary"),
    value = c(epistasis,   Angle,   mutation_rate0,   s,   zeta0,   abserr,   relerr,   err_order,   step_order,   dt0,   itstep,   itermax,   threshold,   unstable_ratio,   unstable_threshold,   rho,   asym,   figures,   family,   density_summary)
  )
  write.csv(params, paste("./", dir.name,"/params.csv",sep=""))
  saveRDS(output, paste("./", dir.name,"/output",sep=""))
  saveRDS(list(f,x.ax,y.ax,t,dt), paste("./", dir.name,"/endstate",sep=""))
} # end play.evol

