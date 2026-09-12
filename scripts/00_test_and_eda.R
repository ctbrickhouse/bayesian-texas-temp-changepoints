# ============================================================================
# 00_test_and_eda.R
#
# Run this first after cloning the repo. It:
#   1. Sanity-checks that the raw data file reads in as expected (a smoke
#      test that your R environment and the data are both set up correctly).
#   2. Cleans the raw data down to summer (Jun-Aug) daily TMAX records and
#      writes data/processed/data_cleaned.csv.
#   3. Runs a short preliminary EDA on the cleaned series.
#
# See README.md > "Running the test / EDA script" for setup instructions.
# ============================================================================

library(here)
library(dplyr)
library(readr)
library(ggplot2)

source(here("R", "clean_data.R"))
source(here("R", "annual_summary.R"))

raw_path   <- here("data", "raw", "dataset_raw.csv")
clean_path <- here("data", "processed", "data_cleaned.csv")

# ---- 1. Smoke tests on the raw file ---------------------------------------

stopifnot("raw data file not found - see README for download instructions" =
            file.exists(raw_path))

raw <- read_csv(raw_path, show_col_types = FALSE)

stopifnot(nrow(raw) > 0)
stopifnot(all(c("STATION", "DATE", "TMAX", "TMAX_ATTRIBUTES") %in% names(raw)))
stopifnot(length(unique(raw$STATION)) == 1)

cat(sprintf(
  "[OK] Raw data: %d rows, %d columns, station %s, dates %s to %s\n",
  nrow(raw), ncol(raw), unique(raw$STATION), min(raw$DATE), max(raw$DATE)
))

# ---- 2. Clean & write data_cleaned.csv -------------------------------------

clean <- clean_summer_tmax(raw_path, out_path = clean_path)

stopifnot(nrow(clean) > 0)
stopifnot(all(clean$month %in% 6:8))
stopifnot(!anyNA(clean$tmax_c))
stopifnot(all(clean$tmax_c > 10 & clean$tmax_c < 50)) # plausible summer TMAX range (C)

cat(sprintf(
  "[OK] Cleaned data written to %s (%d rows, %d-%d)\n",
  clean_path, nrow(clean), min(clean$year), max(clean$year)
))

# ---- 3. Preliminary EDA -----------------------------------------------------

annual <- compute_annual_summer_tmax(clean)

cat("\n-- Data completeness by year --\n")
print(summary(annual$pct_complete))
cat(sprintf(
  "Years meeting >=90%% summer-day completeness: %d of %d (span %d-%d)\n",
  sum(annual$included), nrow(annual), min(annual$year), max(annual$year)
))
if (any(!annual$included)) {
  cat("Excluded (incomplete) years:\n")
  print(annual |> filter(!included) |> select(year, n_valid_days, pct_complete))
} else {
  cat("No years excluded - every year in range has a complete summer record.\n")
}

cat("\n-- Daily summer TMAX (deg F), all included years --\n")
print(summary(filter(clean, year %in% annual$year[annual$included])$tmax_f))

cat("\n-- Annual mean summer TMAX (deg F), first vs. last 10 complete years --\n")
included_years <- sort(annual$year[annual$included])
early <- head(included_years, 10)
late  <- tail(included_years, 10)
print(annual |> filter(year %in% early) |> summarise(mean_tmax_f = mean(mean_tmax_f)))
print(annual |> filter(year %in% late)  |> summarise(mean_tmax_f = mean(mean_tmax_f)))

# 3a. Trend plot: annual mean summer TMAX over time
p1 <- annual |>
  filter(included) |>
  ggplot(aes(year, mean_tmax_f)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE) +
  labs(
    title = "Austin Camp Mabry: Mean Summer (Jun-Aug) Daily Max Temp by Year",
    x = "Year", y = "Mean summer TMAX (deg F)"
  )
print(p1)

# 3b. Distribution shift: earliest vs. most recent complete summers
clean_compare <- clean |>
  mutate(period = case_when(
    year %in% early ~ "Earliest 10 complete yrs",
    year %in% late  ~ "Most recent 10 complete yrs",
    TRUE ~ NA_character_
  )) |>
  filter(!is.na(period))

p2 <- ggplot(clean_compare, aes(tmax_f, fill = period)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Daily Summer TMAX Distribution: Earliest vs. Most Recent Complete Years",
    x = "TMAX (deg F)", fill = NULL
  )
print(p2)

cat("\nAll checks passed and preliminary EDA complete.\n")
