# Bayesian Change-Point Analysis of Texas Summer Temperatures

## Overview
Group project for TAMU STAT 638 Applied Bayesian Analysis (Project 9). See
[`instructions.md`](instructions.md) for the full assignment spec (report
sections, word limits, presentation format).

## Team
- [Campbell Brickhouse](mailto:ctbrickhouse@tamu.edu)
- [Shandon Bringhurst](mailto:shandonb@tamu.edu)
- [Caleb Notheis](mailto:Caleb.Notheis@tamu.edu)

## Data and source
- Daily temperature data from the NOAA National Centers for Environmental
  Information (NCEI), Global Historical Climatology Network-Daily (GHCN-D),
  for the Austin Camp Mabry, Texas station (GHCND: USW00013958).
- The station record begins in 1938 and includes daily maximum temperature
  measurements. Use the longest continuous period with sufficiently complete
  summer observations.
- [Data webpage: NOAA NCEI Daily Summaries — Austin Camp Mabry, TX](https://www.ncei.noaa.gov/cdo-web/datasets/GHCND/stations/GHCND%3AUSW00013958/detail)

## Outcome construction
For each year t, define the summer period as June 1 through August 31. Let
Tmax(t,d) denote the daily maximum temperature on summer day d of year t,
and define Y(t) as the average of the available daily maximum temperatures
over that summer. A year is included only if at least 90% of the summer
days (i.e., at least 83 of 92 days) have a valid, QC-passed daily maximum
temperature; this rule is implemented in `compute_annual_summer_tmax()`
(`R/annual_summary.R`). Preliminary EDA (see below) shows the full
1938-2026 record at this station is in fact complete under this rule, so no
years currently need to be dropped — re-check this if the raw data file is
refreshed.

## Scientific questions
- Has the average summer daily maximum temperature at Austin Camp Mabry
  changed over time?
- Is the long-term pattern better characterized by a gradual trend or by
  one or more abrupt changes? If an abrupt change is supported, when did it
  most likely occur, how large was the change, and how uncertain are the
  estimated timing and magnitude?
- Do conclusions about a change point remain stable under reasonable
  alternative modeling assumptions?

## Repository structure

```
bayesian-texas-temp-changepoints/
├── DESCRIPTION              # package metadata + dependency list (R package management)
├── NAMESPACE                # exports everything in R/
├── .Rprofile                # auto-activates renv if renv/ has been initialized
├── data/
│   ├── raw/
│   │   └── dataset_raw.csv       # raw NOAA GHCN-D download, all variables, all years
│   └── processed/
│       ├── data_cleaned.csv        # cleaned Jun-Aug daily TMAX (from 00_test_and_eda.R)
│       └── data_cleaned_final.csv  # final cleaned data used by the analysis (from 00_data_cleaning_final.R)
├── R/
│   ├── clean_data.R         # clean_summer_tmax(): raw -> cleaned daily TMAX
│   ├── annual_summary.R     # compute_annual_summer_tmax(): daily -> Y(t) with completeness rule
│   └── changepoint_tools.R  # Bayesian models B0-B5: mcp formulas/priors, exact grid posterior of the change point
├── scripts/
│   ├── 00_test_and_eda.R          # smoke tests + builds data_cleaned.csv + EDA
│   ├── 00_data_cleaning_final.R   # step-by-step cleaning -> data_cleaned_final.csv
│   ├── 01_conventional_analysis.R # frequentist analysis (GLS, sup-F, bootstrap) -> analysis/results/
│   ├── 02_bayesian_analysis.R     # Bayesian analysis with mcp -> analysis/results/
│   └── 03_bayesian_sensitivity.R  # prior/AR(1) sensitivity + posterior predictive checks
├── tests/
│   ├── testthat.R           # test runner entry point
│   └── testthat/
│       └── test-clean_data.R
├── analysis/
│   ├── report.Rmd           # written report (reads analysis/results/, cites references.bib)
│   ├── references.bib       # bibliography for the report
│   ├── results/             # saved results and figures from scripts 01-03
│   └── presentation.Rmd     # final slide deck for the live presentation
├── instructions.md          # assignment instructions, unmodified
└── README.md
```

Guidelines for adding to this structure:
- Reusable logic (used by more than one script/Rmd) goes in `R/` as a
  documented function, not copy-pasted.
- One-off / entry-point scripts go in `scripts/`, numbered by run order.
- Anything that reads `data/processed/data_cleaned.csv` should call
  functions from `R/`, not re-implement cleaning logic.
- Model fitting code for the main analysis can live in `R/` (as functions)
  or directly in `analysis/report.Rmd`, whichever keeps the Rmd readable;
  prefer moving anything long/slow into `R/` and calling it from the Rmd.

## Package management

This project is organized like an R package (`DESCRIPTION` + `NAMESPACE` +
`R/`) so dependencies are declared in one place, but it is not built or
installed — scripts just `source()` the files in `R/` directly (see
`scripts/00_test_and_eda.R` for the pattern).

**Quick start (no lockfile yet):**
```r
install.packages("remotes")
remotes::install_deps(dependencies = TRUE)  # installs everything in DESCRIPTION
```

**Recommended: pin exact versions with renv**, so everyone (and the
grader) runs identical package versions:
```r
install.packages("renv")
renv::init()      # first group member only: scans the project, writes renv.lock
```
Commit the resulting `renv.lock` and `renv/` folder (minus `renv/library/`,
already in `.gitignore`) after running this once. Everyone after that just
runs, on each pull:
```r
renv::restore()   # installs the exact pinned versions from renv.lock
```
Any new package used in `R/`, `scripts/`, `tests/`, or `analysis/` should be
added to `DESCRIPTION` (`Imports` for code that must run, `Suggests` for
things only used in tests/reports), then re-run `renv::snapshot()` and
commit the updated `renv.lock`.

## Running the test / EDA script

`scripts/00_test_and_eda.R` is the first thing to run after cloning: it
sanity-checks that the raw data reads in correctly, builds
`data/processed/data_cleaned.csv` from `data/raw/dataset_raw.csv`, and
prints/plots a short preliminary EDA (completeness by year, summary
statistics, an annual-mean trend plot, and an early-vs-recent-years
distribution comparison).

From the project root:
```r
# open the project in RStudio (or set the working directory to the repo root), then:
source("scripts/00_test_and_eda.R")
```
or from a terminal:
```
Rscript scripts/00_test_and_eda.R
```
It will stop with an assertion error if the raw file is missing/malformed
or the cleaned data fails a basic sanity check (non-summer months present,
missing TMAX values, implausible temperatures) — if it completes, your
environment and the data pipeline are both working.

## Running the unit tests

`tests/testthat/test-clean_data.R` contains automated regression tests for
the cleaning and annual-summary functions (column names, no missing
values, unit conversion consistency, completeness-flag logic). Run them
with:
```r
testthat::test_dir("tests/testthat")
```
or, to also re-source `R/` first exactly as the CI/grader would:
```r
source("tests/testthat.R")
```

## Running the analysis

The Bayesian models are fitted with the `mcp` package, which needs **JAGS
4.3.x** installed on your computer, outside R. Get it from
https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/ (Windows: run
`JAGS-4.3.2.exe` with the default settings). JAGS 5.x does **not** work with the
current `rjags` package from CRAN. Check the install with
`library(rjags)`, which should print `Linked to JAGS 4.3.x`.

Run the scripts in order from the project root. Each one saves its results to
`analysis/results/` for the next one and for the report:
```
Rscript scripts/00_data_cleaning_final.R    # ~seconds
Rscript scripts/01_conventional_analysis.R  # ~minutes (bootstrap)
Rscript scripts/02_bayesian_analysis.R      # ~10 minutes (MCMC)
Rscript scripts/03_bayesian_sensitivity.R   # ~10 minutes (MCMC)
```

## Rendering the report and presentation

```r
rmarkdown::render("analysis/report.Rmd")        # written report -> analysis/report.html
rmarkdown::render("analysis/presentation.Rmd")  # slides -> analysis/presentation.html
```
The report reads the saved results in `analysis/results/`, so run scripts
00-03 first (see "Running the analysis"). The presentation (not yet updated)
still reads `data/processed/data_cleaned.csv` from `scripts/00_test_and_eda.R`.

## Modeling task
Develop and justify an appropriate Bayesian time-series/change-point
analysis to address the scientific questions. Compare and assess
reasonable competing descriptions of the temperature series, including
models with no systematic change, a gradual temporal trend, and an abrupt
change when appropriate. Quantify uncertainty in any estimated change
point and in the magnitude of temperature change. Compare the Bayesian
analysis with an appropriate frequentist or alternative
change-point/trend analysis.
