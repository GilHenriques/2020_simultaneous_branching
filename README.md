# Code for simultaneous branching project

This repository contains the code to generate data and to generate the figures for the simultaneous branching project.

## Code to generate data

Code to generate simulations can be found in the folder `code_for_submission`:
- `IBS_code.cpp`: generates IBS
- `ibs_random.cpp`: generates IBS with a random game (for the section on asymmetric games)
- `OSS_code.R`: generates OSS
- `PDE_code.R`: code for PDE model; to run PDE use `run_PDE_code.R` (same folder).

## Analysing data and generating figures

The other folders contain R code to analyse data and create figures.

Raw data files generated from OSS simulations are included as `RDS` files.
Raw data files generated frm IBS and PDE are not included because they are too large, but they can be generated by running the code in `code_for_submission`. However, the repository includes the IBS summary data generated from the raw data (as well as the `R` files used to generate summary data), in the form of `RDS` files. 
All figures can be generated using the included `R` code and the summary data (in the IBS case) or raw data (in the OSS/PDE case).

The figure of the fitness landscape was generated in *Mathematica* using the file `fitness_landscape.nb`.

## More

The derivations for Appendix F (stability to large mutations) are in the *Mathematica* notebook `AppendixF.nb`

Feel free to email me at gilhenris-at-gmail-dot-com for any questions.
