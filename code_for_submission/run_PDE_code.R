source("PDE_code.R")

c_x <- 0.9; d_x = 1.65; b_x = 1.05;
b_y = b_x; c_y = c_x; d_y = d_x;

# single replicate run: set to TRUE
if (TRUE) {
  start_time <- Sys.time()
  NAME <- "test_run_225"
  if(sum(list.files()==NAME) == 0){

    play.evol(
      epistasis = "additive",
      Angle = 225,
      name = NAME,
      itstep = 5,
      mutation_rate0 = 0.001,
      zeta0 = 2e-7,
      s = 100,
      itermax = 2500,
      family = "quadratic",
      density_summary = FALSE,
      figures = TRUE,
      asym = 0.5    )
  }
  end_time <- Sys.time()
  end_time - start_time
}
