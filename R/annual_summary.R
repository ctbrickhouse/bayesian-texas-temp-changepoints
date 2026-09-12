#' Construct the annual summer TMAX outcome Y(t)
#'
#' For each year, averages the available daily TMAX values over June 1 -
#' August 31 and flags whether the year meets a data-completeness threshold
#' (default: at least 90% of the 92 summer days present with valid TMAX,
#' per the project's outcome-construction rule).
#'
#' @param clean_df Cleaned daily data frame, as returned by
#'   \code{\link{clean_summer_tmax}}.
#' @param days_per_summer Number of calendar days in Jun 1 - Aug 31
#'   (30 + 31 + 31 = 92; this does not vary with leap years since Feb is
#'   not in the summer window).
#' @param completeness_threshold Minimum fraction of \code{days_per_summer}
#'   that must have a valid TMAX for a year to be included in Y(t).
#'
#' @return A tibble with one row per year: year, n_valid_days,
#'   pct_complete, included, mean_tmax_c, mean_tmax_f.
compute_annual_summer_tmax <- function(clean_df,
                                        days_per_summer = 92,
                                        completeness_threshold = 0.90) {
  clean_df |>
    dplyr::group_by(year) |>
    dplyr::summarise(
      n_valid_days = dplyr::n(),
      mean_tmax_c = mean(tmax_c),
      mean_tmax_f = mean(tmax_f),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      pct_complete = n_valid_days / days_per_summer,
      included = pct_complete >= completeness_threshold
    ) |>
    dplyr::arrange(year) |>
    dplyr::select(year, n_valid_days, pct_complete, included, mean_tmax_c, mean_tmax_f)
}
