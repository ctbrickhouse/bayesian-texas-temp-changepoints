# This project is organized like an R package but is not installed or
# R CMD check-ed; tests are run directly against the sourced R/ functions.
library(testthat)

source(here::here("R", "clean_data.R"))
source(here::here("R", "annual_summary.R"))

test_dir(here::here("tests", "testthat"))
