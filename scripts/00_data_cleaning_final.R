library(here)
library(dplyr)
library(readr)
library(lubridate)

# This script does the same thing as clean_summer_tmax() in R/clean_data.R,
# one small step at a time, so you can see exactly what happens to the data.
# Run it line by line and look at the result after each step, e.g.
# with View(step1) or head(step1).

# ---- 0. load the raw data -------------------------------------
raw_path   <- here("data", "raw", "dataset_raw.csv")
raw <- read_csv(raw_path, show_col_types = FALSE)

nrow(raw)   # number of days in the raw file (all months, all years)

# ---- 1. clean the raw data -------------------------------------

# select relevant columns
# The raw file has ~150 columns; we only need these four:
#   STATION         = weather station ID (Austin Camp Mabry)
#   DATE            = the day of the observation
#   TMAX            = maximum temperature that day, in TENTHS of a degree C
#                     (e.g. 328 means 32.8 C)
#   TMAX_ATTRIBUTES = flags describing the TMAX value, e.g. ",,0"
step1 <- raw |>
  select(STATION, DATE, TMAX, TMAX_ATTRIBUTES)

# split the date into year, month and day columns
step2 <- step1 |>
  mutate(
    year  = year(DATE),
    month = month(DATE),
    day   = day(DATE)
  )

# keep only summer days (June = 6, July = 7, August = 8)
step3 <- step2 |>
  filter(month %in% c(6, 7, 8))

nrow(step3)   # number of summer days

# drop days where the temperature is missing (there aren't any)
step4 <- step3 |>
  filter(!is.na(TMAX))

nrow(step4)   # summer days that have a temperature

# pull out the quality-control (QC) flag
# TMAX_ATTRIBUTES looks like "M,Q,S": measurement flag, QC flag, source flag.
# The QC flag is the 2nd piece. It is empty ("") when the value passed
# NOAA's quality checks, and a letter when the value failed a check.
step5 <- step4 |>
  mutate(
    mflag = sapply(strsplit(TMAX_ATTRIBUTES, ","), function(x) x[1]),
    qflag = sapply(strsplit(TMAX_ATTRIBUTES, ","), function(x) x[2]),
    sflag = sapply(strsplit(TMAX_ATTRIBUTES, ","), function(x) x[3])
  )

# # keep only days that passed quality control and have no measurement flags (no rows are removed)
step6 <- step5 |>
  filter(qflag == "" & mflag == "")

nrow(step6)

# convert temperature to degrees Celsius
step7 <- step6 |>
  mutate(
    tmax_c = round(TMAX / 10, 2)
  )

# keep the final columns, give them lowercase names, and sort by date
clean <- step7 |>
  select(
    station = STATION,
    date    = DATE,
    year, month, day,
    tmax_c, sflag
  ) |>
  arrange(date)

# ---- 2. check and save the cleaned data ------------------------

head(clean)
summary(clean$tmax_c)
range(clean$year)

write_csv(clean, here("data", "processed", "data_cleaned_final.csv"))
