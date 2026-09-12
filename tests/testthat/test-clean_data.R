test_that("clean_summer_tmax filters to summer months and drops missing/QC-failed values", {
  raw_path <- here::here("data", "raw", "dataset_raw.csv")
  skip_if_not(file.exists(raw_path), "raw data file not present")

  clean <- clean_summer_tmax(raw_path)

  expect_true(nrow(clean) > 0)
  expect_true(all(clean$month %in% 6:8))
  expect_false(anyNA(clean$tmax_c))
  expect_false(anyNA(clean$tmax_f))
  expect_true(all(clean$qc_flag == ""))
  expect_setequal(
    names(clean),
    c("station", "date", "year", "month", "day", "tmax_c", "tmax_f", "qc_flag")
  )
  # tmax_f should be a consistent linear transform of tmax_c
  expect_equal(clean$tmax_f, round(clean$tmax_c * 9 / 5 + 32, 2))
})

test_that("clean_summer_tmax writes a CSV when out_path is supplied", {
  raw_path <- here::here("data", "raw", "dataset_raw.csv")
  skip_if_not(file.exists(raw_path), "raw data file not present")

  tmp <- withr::local_tempfile(fileext = ".csv")
  clean <- clean_summer_tmax(raw_path, out_path = tmp)

  expect_true(file.exists(tmp))
  reread <- readr::read_csv(tmp, show_col_types = FALSE)
  expect_equal(nrow(reread), nrow(clean))
})

test_that("compute_annual_summer_tmax flags completeness correctly", {
  fake <- data.frame(
    year = c(rep(2020, 92), rep(2021, 40)),
    tmax_c = 30,
    tmax_f = 86
  )

  annual <- compute_annual_summer_tmax(fake)

  expect_equal(nrow(annual), 2)
  expect_true(annual$included[annual$year == 2020])
  expect_false(annual$included[annual$year == 2021])
  expect_equal(annual$pct_complete[annual$year == 2021], 40 / 92)
})
