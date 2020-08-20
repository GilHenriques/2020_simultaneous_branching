library(tidyverse)

# This is a ridiculously slow way to calculate a regression line,
# but I only need to run it once.
# I use it because the axes are circular, i.e. angle = pi is the same as angle = 0.
# So a normal regression doesn't work.

angles <- readRDS("2020-07 IBS/angles_225_frequent_output_rds")

output <- rbindlist(angles, use.names = TRUE) %>% as_tibble()

angle_onset_branching <- output %>% 
  filter(end_branching | onset) %>% 
  filter(nr_clusters == 2)

angle_before_onset <- filter(output, nr_clusters == 2)[which(filter(output, nr_clusters == 2)$onset)-1, ]$angle
angle_branching <- angle_onset_branching$angle[angle_onset_branching$end_branching]


plot(angle_branching ~ angle_before_onset)

residual_squared <- function(x, y, a, b){
  yhat <- a*x + b
  square1 <- (yhat - y)^2
  square2 <- (yhat + pi - y)^2
  square3 <- (yhat - pi - y)^2
  return(min(square1, square2, square3))
}

a <- 1
b <- 0
sd <- 0.0001
SS_old <- 10000
output <- tibble(counter = 99, SS = 99, a = 99, b = 99)
iterations <- 50000
for(counter in 1:iterations){
  SS <- 0
  a_mut <- rnorm(1, a, sd)
  b_mut <- rnorm(1, b, sd)
  for(i in 1:length(angle_before_onset)){
    SS <- SS + residual_squared(angle_before_onset[i], angle_branching[i], a_mut, b_mut)
  }
  if(SS < SS_old) {
    a <- a_mut
    b <- b_mut
    output <- rbind(output, c(counter, SS, a, b))
    SS_old <- SS
  }
  counter <- counter+1
  if(counter %% 100 == 0) print(paste0(round(counter/iterations*100,3), " % complete"))
}
output <- output[-1,]
par(mfrow = c(1, 3))
plot(output$SS ~ output$counter)
plot(output$a ~ output$counter)
plot(output$b ~ output$counter)
par(mfrow = c(1,1))
output %>% tail

plot(angle_branching ~ angle_before_onset)
lines(seq(from = 0, to = pi, length = 1000), output[nrow(output),]$a * seq(from = 0, to = pi, length = 1000) + output[nrow(output),]$b, col = "red", lwd = 5)
output[nrow(output),]$a # 1.029445
output[nrow(output),]$b # -0.1378501


n <- length(angle_before_onset)
df <- n - 2
tstar <- abs(qt(0.05/2, df))

sy <- sqrt(output[nrow(output),]$SS/(n-2))
Xvec <- seq(from = 0, to = pi, length = 1000)
Xbar <- mean(angle_before_onset)

A <- output[nrow(output),]$a 
B <-  output[nrow(output),]$b
Yhat <-  A*Xvec + B
conf_top <- Yhat + tstar * sy * sqrt((1/n)+ ((Xvec-Xbar)^2)/sum((Xvec-Xbar)^2))
conf_bot <- Yhat - tstar * sy * sqrt((1/n)+ ((Xvec-Xbar)^2)/sum((Xvec-Xbar)^2))

plot(angle_branching ~ angle_before_onset)
lines(Xvec, Yhat, col = "red", lwd = 3)
lines(Xvec, conf_top, col = "red", lwd = 2)
lines(Xvec, conf_bot, col = "red", lwd = 2)

saveRDS(object = tibble(x = Xvec, y = Yhat, ymin = conf_bot, ymax = conf_top), file = "2020-07 IBS/regression_line")
