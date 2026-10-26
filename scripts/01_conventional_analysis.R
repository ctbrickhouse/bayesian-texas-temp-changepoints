library(here)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(nlme)
library(strucchange)
library(trend)

# Conventional (frequentist) analysis of Y(t) = average summer (Jun-Aug)
# daily maximum temperature at Austin Camp Mabry, 1938-2026.
#
# Candidate mean models (s = year - 1982, the middle year):
#   M0  y = a0                                      (no change)
#   M1  y = a0 + a1*s                               (linear trend)
#   M2  y = a0 + a1*s + a2*s^2                      (quadratic trend)
#   M3  y = a0 + a1*s + b0*I(year > g)              (jump in level)
#   M4  y = a0 + a1*s + b1*(year - g)+              (change in slope, continuous)
#   M5  y = a0 + a1*s + b0*I(year > g) + b1*(year - g)+   (jump + slope change)
# Each one is fit with iid errors and with AR(1) errors
#   e(t) = rho*e(t-1) + u(t),  u(t) ~ N(0, sigma_u^2).
# All fits use maximum likelihood (not REML) so that log-likelihoods are
# comparable across different mean models.

set.seed(638)
B <- 300   # bootstrap replicates (raise for the final report)

out_dir <- here("analysis", "results")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

col_data <- "grey55"
col_m1   <- "#2a78d6"
col_cp   <- "#eb6834"

# ---- 0. build Y(t) ---------------------------------------------

daily <- read_csv(here("data", "processed", "data_cleaned_final.csv"),
                  show_col_types = FALSE)

annual <- daily |>
  group_by(year) |>
  summarise(n_days = n(), y = mean(tmax_c), .groups = "drop") |>
  arrange(year)

# every summer has all 92 days, so no year needs to be dropped
stopifnot(all(annual$n_days == 92), all(diff(annual$year) == 1))

mid_year <- (min(annual$year) + max(annual$year)) / 2   # 1982
annual <- annual |> mutate(s = year - mid_year)
n <- nrow(annual)

# candidate change-point years: trim 15% of the record at each end
gamma_grid <- seq(min(annual$year) + ceiling(0.15 * n),
                  max(annual$year) - ceiling(0.15 * n))

# ---- 1. helpers ------------------------------------------------

cp_formula <- list(
  M3 = y ~ s + step,
  M4 = y ~ s + hinge,
  M5 = y ~ s + step + hinge
)

add_cp_terms <- function(d, g) {
  d |> mutate(step = as.numeric(year > g), hinge = pmax(year - g, 0))
}

fit_gls <- function(formula, d, ar1) {
  gls(formula, data = d, method = "ML",
      correlation = if (ar1) corAR1(form = ~ year) else NULL)
}

# fit a change-point model at a fixed change year g
fit_cp <- function(g, model, d, ar1) {
  fit_gls(cp_formula[[model]], add_cp_terms(d, g), ar1)
}

# profile the log-likelihood over the grid and keep the best change year
profile_cp <- function(model, d, ar1, grid = gamma_grid) {
  ll <- vapply(grid, function(g) {
    f <- tryCatch(fit_cp(g, model, d, ar1), error = function(e) NULL)
    if (is.null(f)) -Inf else as.numeric(logLik(f))
  }, numeric(1))
  g_hat <- grid[which.max(ll)]
  list(gamma = g_hat, profile = tibble(gamma = grid, loglik = ll),
       fit = fit_cp(g_hat, model, d, ar1))
}

# information criteria; change-point models get one extra parameter for gamma
ic_row <- function(fit, model, ar1, has_gamma) {
  ll <- as.numeric(logLik(fit))
  k  <- attr(logLik(fit), "df") + has_gamma
  m  <- fit$dims$N
  tibble(model = model, errors = if (ar1) "AR(1)" else "iid",
         logLik = ll, k = k,
         AIC = -2 * ll + 2 * k,
         AICc = -2 * ll + 2 * k + 2 * k * (k + 1) / (m - k - 1),
         BIC = -2 * ll + log(m) * k)
}

get_rho <- function(fit) {
  if (is.null(fit$modelStruct$corStruct)) return(0)
  unname(coef(fit$modelStruct$corStruct, unconstrained = FALSE))
}

# ---- 2. EDA ----------------------------------------------------

ols <- lm(y ~ s, data = annual)
acf_ols <- acf(resid(ols), plot = FALSE, lag.max = 15)
cat("Lag-1 autocorrelation of OLS residuals:", round(acf_ols$acf[2], 3), "\n")

# ---- 3. fit all 12 models --------------------------------------

poly_formula <- list(M0 = y ~ 1, M1 = y ~ s, M2 = y ~ s + I(s^2))

fits <- list()
profiles <- list()
for (ar1 in c(FALSE, TRUE)) {
  tag <- if (ar1) "AR1" else "iid"
  for (m in names(poly_formula)) {
    fits[[paste(m, tag)]] <- fit_gls(poly_formula[[m]], annual, ar1)
  }
  for (m in names(cp_formula)) {
    p <- profile_cp(m, annual, ar1)
    fits[[paste(m, tag)]] <- p$fit
    profiles[[paste(m, tag)]] <- p$profile |> mutate(model = m, errors = tag)
    attr(fits[[paste(m, tag)]], "gamma") <- p$gamma
  }
}

# ---- 4. model comparison ---------------------------------------

comparison <- bind_rows(lapply(names(fits), function(nm) {
  f <- fits[[nm]]
  parts <- strsplit(nm, " ")[[1]]
  ic_row(f, parts[1], parts[2] == "AR1", has_gamma = parts[1] %in% names(cp_formula)) |>
    mutate(gamma = if (is.null(attr(f, "gamma"))) NA_real_ else attr(f, "gamma"),
           rho = get_rho(f),
           sigma = f$sigma)
})) |>
  mutate(dAICc = AICc - min(AICc), dBIC = BIC - min(BIC)) |>
  arrange(AICc)

print(comparison, width = Inf)

# likelihood-ratio tests for nested models (gamma fixed models excluded)
lrt <- function(small, big) {
  stat <- 2 * (as.numeric(logLik(fits[[big]])) - as.numeric(logLik(fits[[small]])))
  df <- attr(logLik(fits[[big]]), "df") - attr(logLik(fits[[small]]), "df")
  tibble(test = paste(small, "vs", big), LR = stat, df = df,
         p = pchisq(stat, df, lower.tail = FALSE))
}
lrts <- bind_rows(
  lrt("M0 AR1", "M1 AR1"), lrt("M1 AR1", "M2 AR1"),
  lrt("M0 iid", "M0 AR1"), lrt("M1 iid", "M1 AR1"), lrt("M2 iid", "M2 AR1")
)
print(lrts)

# Change-point tests. gamma does not exist under "no change", so the usual
# chi-square LRT is invalid. Two alternatives:
# (a) Andrews sup-F test (assumes iid errors)
fs <- Fstats(y ~ s, data = annual, from = 0.15)
supF <- sctest(fs, type = "supF")
print(supF)

# (b) parametric bootstrap of the sup-LR statistic, simulating from the
#     fitted M1 + AR(1) null model so the autocorrelation is respected
null_fit <- fits[["M1 AR1"]]
rho0 <- get_rho(null_fit)
sd_u0 <- null_fit$sigma * sqrt(1 - rho0^2)   # innovation sd
mu0 <- fitted(null_fit)

sup_lr <- function(d, model) {
  ll0 <- as.numeric(logLik(fit_gls(y ~ s, d, TRUE)))
  ll1 <- max(profile_cp(model, d, TRUE)$profile$loglik)
  2 * (ll1 - ll0)
}

cp_test_models <- c("M3", "M4", "M5")
obs_sup_lr <- sapply(cp_test_models, function(m) sup_lr(annual, m))

boot_sup_lr <- replicate(B, {
  d <- annual
  d$y <- mu0 + as.numeric(arima.sim(list(ar = rho0), n = n, sd = sd_u0))
  sapply(cp_test_models, function(m) tryCatch(sup_lr(d, m), error = function(e) NA))
})
cp_tests <- tibble(
  model = cp_test_models,
  sup_LR = obs_sup_lr,
  boot_p = sapply(cp_test_models, function(m) mean(boot_sup_lr[m, ] >= obs_sup_lr[m], na.rm = TRUE))
)
print(cp_tests)

# ---- 5. parameter estimates for the leading models -------------

best_cp <- comparison |> filter(model %in% names(cp_formula)) |> slice(1)
best_cp_name <- paste(best_cp$model, if (best_cp$errors == "AR(1)") "AR1" else "iid")
key_models <- unique(c("M1 iid", "M1 AR1", "M2 AR1", best_cp_name))

estimates <- bind_rows(lapply(key_models, function(nm) {
  f <- fits[[nm]]
  ci <- intervals(f, which = "coef")$coef
  tibble(fit = nm, term = rownames(ci), estimate = ci[, "est."],
         lower = ci[, "lower"], upper = ci[, "upper"],
         se = sqrt(diag(vcov(f))))
}))
print(estimates, n = Inf)

# warming rate in degrees C per decade from the linear-trend models
trend_rate <- estimates |>
  filter(fit %in% c("M1 iid", "M1 AR1"), term == "s") |>
  mutate(across(c(estimate, lower, upper, se), ~ .x * 10))
print(trend_rate)

rho_ci <- intervals(fits[["M1 AR1"]], which = "var-cov")$corStruct
rho_hat <- rho_ci[, "est."]
n_eff <- n * (1 - rho_hat) / (1 + rho_hat)
cat("rho (M1 AR1):", round(rho_ci, 3), " effective n ~", round(n_eff, 1), "\n")

# ---- 6. residual diagnostics -----------------------------------

diag_tbl <- bind_rows(lapply(key_models, function(nm) {
  r <- resid(fits[[nm]], type = "normalized")
  fitdf <- as.integer(grepl("AR1", nm))
  tibble(fit = nm,
         ljung_box_p_lag5  = Box.test(r, lag = 5,  type = "Ljung-Box", fitdf = fitdf)$p.value,
         ljung_box_p_lag10 = Box.test(r, lag = 10, type = "Ljung-Box", fitdf = fitdf)$p.value,
         shapiro_p = shapiro.test(r)$p.value)
}))
print(diag_tbl)

# Sensitivity: richer error processes on the linear-trend and best change-point models
arma_sens <- bind_rows(lapply(c("M1", best_cp$model), function(m) {
  d <- if (m == "M1") annual else add_cp_terms(annual, best_cp$gamma)
  fmla <- if (m == "M1") y ~ s else cp_formula[[m]]
  errs <- list("AR(1)"     = corAR1(form = ~ year),
               "AR(2)"     = corARMA(form = ~ year, p = 2, q = 0),
               "ARMA(1,1)" = corARMA(form = ~ year, p = 1, q = 1))
  bind_rows(lapply(names(errs), function(e) {
    f <- gls(fmla, data = d, method = "ML", correlation = errs[[e]])
    tibble(model = m, errors = e, AIC = AIC(f),
           slope_per_decade = 10 * coef(f)[["s"]],
           slope_se_decade = 10 * sqrt(vcov(f)["s", "s"]))
  }))
}))
print(arma_sens)

# ---- 7. change-point uncertainty (parametric bootstrap) --------

cp_fit <- fits[[best_cp_name]]
cp_ar1 <- best_cp$errors == "AR(1)"
rho_cp <- get_rho(cp_fit)
sd_u_cp <- cp_fit$sigma * sqrt(1 - rho_cp^2)
mu_cp <- fitted(cp_fit)

boot_gamma <- replicate(B, {
  d <- annual
  d$y <- mu_cp + as.numeric(
    if (rho_cp > 0) arima.sim(list(ar = rho_cp), n = n, sd = sd_u_cp) else rnorm(n, 0, sd_u_cp)
  )
  profile_cp(best_cp$model, d, cp_ar1)$gamma
})
gamma_ci <- quantile(boot_gamma, c(0.025, 0.5, 0.975))
cat("Best change-point model:", best_cp_name, " gamma-hat =", best_cp$gamma,
    " bootstrap 95% interval:", gamma_ci[c(1, 3)], "\n")

# ---- 8. leave-future-out check (train 1938-2006, test 2007-2026) -----

split_year <- 2006
train <- annual |> filter(year <= split_year)
test  <- annual |> filter(year >  split_year)
train_grid <- gamma_grid[gamma_grid <= split_year - ceiling(0.15 * nrow(train))]

forecast_model <- function(m, ar1) {
  if (m %in% names(cp_formula)) {
    p <- profile_cp(m, train, ar1, grid = train_grid)
    f <- p$fit
    newd <- add_cp_terms(test, p$gamma)
  } else {
    f <- fit_gls(poly_formula[[m]], train, ar1)
    newd <- test
  }
  mu <- as.numeric(predict(f, newdata = newd))
  # AR(1): the last training residual decays as rho^h
  e_last <- tail(resid(f), 1)
  h <- test$year - split_year
  pred <- mu + get_rho(f)^h * e_last
  tibble(model = m, errors = if (ar1) "AR(1)" else "iid",
         RMSE = sqrt(mean((test$y - pred)^2)),
         MAE  = mean(abs(test$y - pred)))
}
lfo <- bind_rows(lapply(c(FALSE, TRUE), function(ar1)
  bind_rows(lapply(c(names(poly_formula), names(cp_formula)), forecast_model, ar1 = ar1)))) |>
  arrange(RMSE)
print(lfo)

# ---- 9. nonparametric cross-check: Mann-Kendall and Sen slope -----

mk_raw <- mk.test(annual$y)
sen_raw <- sens.slope(annual$y)

# prewhitening (von Storch): remove lag-1 autocorrelation of the detrended
# series before testing, so the test is not fooled by persistence
r1 <- acf(resid(ols), plot = FALSE)$acf[2]
y_pw <- annual$y[-1] - r1 * annual$y[-n]
mk_pw <- mk.test(y_pw)
sen_pw <- sens.slope(annual$y[-1] - r1 * annual$y[-n])

nonparam <- tibble(
  version = c("raw", "prewhitened"),
  mk_p = c(mk_raw$p.value, mk_pw$p.value),
  sen_slope_per_decade = 10 * c(sen_raw$estimates, sen_pw$estimates / (1 - r1)),
  sen_lower_decade = 10 * c(sen_raw$conf.int[1], sen_pw$conf.int[1] / (1 - r1)),
  sen_upper_decade = 10 * c(sen_raw$conf.int[2], sen_pw$conf.int[2] / (1 - r1))
)
print(nonparam)

# ---- 10. figures -----------------------------------------------

fit_lines <- bind_rows(
  tibble(year = annual$year, value = fitted(fits[["M1 AR1"]]), fit = "M1 linear trend, AR(1)"),
  tibble(year = annual$year, value = fitted(cp_fit),
         fit = paste0(best_cp$model, " change point (", best_cp$gamma, "), ", best_cp$errors))
)

p_fit <- ggplot(annual, aes(year, y)) +
  geom_point(colour = col_data, size = 1.8) +
  geom_line(data = fit_lines, aes(year, value, colour = fit), linewidth = 0.9) +
  scale_colour_manual(values = c(col_m1, col_cp), name = NULL) +
  labs(x = NULL, y = "Mean summer TMAX (°C)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "conventional_fits.png"), p_fit, width = 8, height = 4.5, dpi = 150)

p_prof <- bind_rows(profiles) |>
  group_by(model, errors) |>
  mutate(rel_loglik = loglik - max(loglik)) |>
  ggplot(aes(gamma, rel_loglik)) +
  geom_line(colour = col_m1, linewidth = 0.8) +
  geom_hline(yintercept = -qchisq(0.95, 1) / 2, linetype = "dashed", colour = "grey50") +
  facet_grid(errors ~ model) +
  labs(x = "Candidate change year", y = "Profile log-likelihood (relative to max)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "conventional_cp_profiles.png"), p_prof, width = 9, height = 5, dpi = 150)

acf_df <- bind_rows(lapply(c("M1 iid", "M1 AR1"), function(nm) {
  a <- acf(resid(fits[[nm]], type = "normalized"), plot = FALSE, lag.max = 15)
  tibble(fit = nm, lag = a$lag[-1], acf = a$acf[-1])
}))
p_acf <- ggplot(acf_df, aes(lag, acf)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_hline(yintercept = c(-1, 1) * 1.96 / sqrt(n), linetype = "dashed", colour = "grey50") +
  geom_segment(aes(xend = lag, yend = 0), colour = col_m1, linewidth = 0.9) +
  facet_wrap(~ fit) +
  labs(x = "Lag (years)", y = "ACF of normalized residuals") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "conventional_resid_acf.png"), p_acf, width = 8, height = 3.5, dpi = 150)

# ---- 11. save results for the report ---------------------------

results <- list(
  annual = annual, comparison = comparison, lrts = lrts, supF = supF,
  cp_tests = cp_tests, estimates = estimates, trend_rate = trend_rate,
  rho_ci = rho_ci, n_eff = n_eff, diagnostics = diag_tbl, arma_sens = arma_sens,
  best_cp = best_cp, boot_gamma = boot_gamma, gamma_ci = gamma_ci,
  lfo = lfo, nonparam = nonparam, profiles = bind_rows(profiles), B = B
)
saveRDS(results, file.path(out_dir, "conventional_results.rds"))
write_csv(comparison, file.path(out_dir, "conventional_model_comparison.csv"))
