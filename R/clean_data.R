#' Clean the raw GHCN-D daily file down to summer (Jun-Aug) TMAX records
#'
#' Reads the raw NOAA GHCN-D daily summaries CSV for Austin Camp Mabry
#' (GHCND:USW00013958), keeps only June-August observations, drops
#' observations with a missing or QC-failed TMAX value, converts TMAX from
#' tenths of degrees Celsius to Celsius and Fahrenheit, and returns a tidy
#' data frame with only the columns needed for the change-point analysis.
#'
#' @param raw_path Path to the raw CSV (as downloaded from NOAA NCEI).
#' @param out_path Optional path to write the cleaned CSV to. If NULL
#'   (default), nothing is written and only the data frame is returned.
#'
#' @return A tibble with one row per summer day and columns: station, date,
#'   year, month, day, tmax_c, tmax_f, qc_flag.
clean_summer_tmax <- function(raw_path, out_path = NULL) {
  raw <- readr::read_csv(
    raw_path,
    col_select = c("STATION", "DATE", "TMAX", "TMAX_ATTRIBUTES"),
    col_types = readr::cols(
      STATION = readr::col_character(),
      DATE = readr::col_date(format = "%Y-%m-%d"),
      TMAX = readr::col_double(),
      TMAX_ATTRIBUTES = readr::col_character()
    )
  )

  cleaned <- raw |>
    dplyr::mutate(
      year = lubridate::year(DATE),
      month = lubridate::month(DATE),
      day = lubridate::day(DATE),
      # TMAX_ATTRIBUTES format: "Measurement flag,Quality flag,Source flag[,Time]"
      qc_flag = vapply(
        strsplit(TMAX_ATTRIBUTES, ",", fixed = TRUE),
        function(x) if (length(x) >= 2) x[2] else "",
        character(1)
      )
    ) |>
    dplyr::filter(
      month %in% 6:8,
      !is.na(TMAX),
      qc_flag == ""
    ) |>
    dplyr::mutate(
      tmax_c = round(TMAX / 10, 2),
      tmax_f = round(tmax_c * 9 / 5 + 32, 2)
    ) |>
    dplyr::transmute(
      station = STATION,
      date = DATE,
      year, month, day,
      tmax_c, tmax_f, qc_flag
    ) |>
    dplyr::arrange(date)

  if (!is.null(out_path)) {
    dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(cleaned, out_path)
  }

  cleaned
}
