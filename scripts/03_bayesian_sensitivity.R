library(here)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(mcp)
library(posterior)
library(bayesplot)

source(here("R", "changepoint_tools.R"))

# Sensitivity analysis and model checking for the Bayesian analysis.
# Run 02_bayesian_analysis.R first: this script re-uses its settings and draws.
#
#   1. Prior predictive check: do the priors produce plausible summers?
#   2. Prior sensitivity of the jump model B3 (change-point window, tighter
#      and wider priors).
#   3. Error-model sensitivity: AR(1) errors for B1 and B3.
#   4. Posterior predictive checks for B1, B3 and B4.

set.seed(6380)

out_dir <- here("analysis", "results")
fig_dir <- file.path(out_dir, "figures")

col_data <- "grey55"
col_cp   <- "#eb6834"

bayes <- readRDS(file.path(out_dir, "bayesian_results.rds"))
conv  <- readRDS(file.path(out_dir, "conventional_results.rds"))
annual <- conv$annual
d <- annual |> select(s, y) |> as.data.frame()

# re-use exactly the same settings as script 02
st <- bayes$settings
mcmc <- st$mcmc
mid_year <- st$mid_year
s_range <- range(d$s)

# ---- 1. prior predictive check ---------------------------------

# Simulate whole temperature series from the B3 prior (no data involved) and
# look at the summer means they imply. Weakly informative priors should allow
# every physically plausible value but rule out absurd ones.
n_sim <- 4000
prior_draws <- tibble(
  cp_1    = runif(n_sim, st$cp_lower, st$cp_upper),
  int_1   = rnorm(n_sim, st$prior_mean[["int_1"]], st$prior_sd[["int_1"]]),
  int_2   = rnorm(n_sim, st$prior_mean[["int_2"]], st$prior_sd[["int_2"]]),
  s_1     = rnorm(n_sim, st$prior_mean[["s_1"]], st$prior_sd[["s_1"]]),
  sigma_1 = abs(rnorm(n_sim, 0, st$sigma_sd))
)
mu_prior <- mean_curve("B3", prior_draws, d$s)
y_prior  <- mu_prior + matrix(rnorm(length(mu_prior), 0, prior_draws$sigma_1), nrow(mu_prior))
prior_predictive <- tibble(
  lower_2.5 = quantile(y_prior, 0.025),
  upper_97.5 = quantile(y_prior, 0.975),
  share_below_25 = mean(y_prior < 25),
  share_above_45 = mean(y_prior > 45),
  observed_min = min(d$y),
  observed_max = max(d$y)
)
print(prior_predictive)

# ---- 2. prior sensitivity of the jump model B3 -----------------

# Each variant changes one thing relative to the main analysis in script 02.
variants <- list(
  "Main analysis" = list(),
  "Change point anywhere in 1939-2025" = list(cp_lower = -43, cp_upper = 43),
  "Priors 4x tighter" = list(prior_sd = st$prior_sd / 4),
  "Priors 4x wider" = list(prior_sd = st$prior_sd * 4),
  "AR(1) errors" = list(ar = TRUE)
)

run_variant <- function(v) {
  args <- modifyList(list(prior_mean = st$prior_mean, prior_sd = st$prior_sd,
                          sigma_sd = st$sigma_sd, cp_lower = st$cp_lower,
                          cp_upper = st$cp_upper, ar = FALSE), v)
  fit_cp_by_region("B3", d, args$prior_mean, args$prior_sd, args$sigma_sd,
                   args$cp_lower, args$cp_upper, split = st$split, mcmc = mcmc,
                   ar = args$ar, rho_grid = seq(-0.9, 0.9, by = 0.05))
}

sens_runs <- list()
for (v in names(variants)) {
  cat("\n==== B3 variant:", v, "====\n")
  if (v == "Main analysis") {
    # already fitted in script 02
    sens_runs[[v]] <- list(draws = bayes$draws$B3, weights = bayes$region_weights$B3)
  } else {
    sens_runs[[v]] <- run_variant(variants[[v]])
    sens_runs[[v]]$draws <- add_derived(sens_runs[[v]]$draws, "B3", s_range, mid_year)
  }
}

# R-hat / ESS of the new MCMC runs
sens_convergence <- bind_rows(lapply(setdiff(names(sens_runs), "Main analysis"), function(v) {
  bind_rows(lapply(c("early", "late"), function(r) {
    summarise_draws(as_draws_df(sens_runs[[v]]$fits[[r]]$mcmc_post), "rhat", "ess_bulk") |>
      filter(variable != "s_2") |>
      mutate(variant = v, region = r, .before = 1)
  }))
}))
cat("Sensitivity runs - largest R-hat:", round(max(sens_convergence$rhat), 4),
    " smallest bulk ESS:", round(min(sens_convergence$ess_bulk)), "\n")

# one row per variant with the quantities the conclusions rest on
sensitivity <- bind_rows(lapply(names(sens_runs), function(v) {
  dr <- sens_runs[[v]]$draws
  late <- filter(dr, region == "late")
  early <- filter(dr, region == "early")
  tw <- summ(dr$total_warming)
  tibble(
    variant = v,
    p_late = sens_runs[[v]]$weights[["late"]],          # P(change after 1987)
    cp_late_median = median(late$cp_year),
    cp_early_median = median(early$cp_year),
    jump_late = mean(late$jump),
    jump_late_lower = unname(quantile(late$jump, 0.025)),
    jump_late_upper = unname(quantile(late$jump, 0.975)),
    jump_early = mean(early$jump),
    total_warming = tw$mean, tw_lower = tw$lower, tw_upper = tw$upper,
    rho = if ("ar1_1" %in% names(dr)) mean(dr$ar1_1) else NA_real_,
    rho_lower = if ("ar1_1" %in% names(dr)) unname(quantile(dr$ar1_1, 0.025)) else NA_real_,
    rho_upper = if ("ar1_1" %in% names(dr)) unname(quantile(dr$ar1_1, 0.975)) else NA_real_
  )
}))
print(sensitivity, width = Inf)

# ---- 3. AR(1) errors for the linear-trend model B1 -------------

# B1 has no change point, so a single mcp run is enough
b1_ar_fit <- fit_mcp("B1", d,
                     mcp_priors("B1", st$prior_mean, st$prior_sd, st$sigma_sd, ar = TRUE),
                     mcmc, ar = TRUE)
b1_ar <- add_derived(get_draws(b1_ar_fit), "B1", s_range, mid_year)
b1_ar_summary <- bind_rows(
  summ(bayes$draws$B1$trend_decade) |> mutate(errors = "iid", quantity = "trend_decade"),
  summ(b1_ar$trend_decade) |> mutate(errors = "AR(1)", quantity = "trend_decade"),
  summ(b1_ar$ar1_1) |> mutate(errors = "AR(1)", quantity = "rho")
) |> relocate(errors, quantity)
print(b1_ar_summary)

# ---- 4. posterior predictive checks ----------------------------

# Replicated data sets y_rep: for each of 4000 posterior draws, simulate a new
# 1938-2026 series from the model. If the model is adequate, the observed
# series should look like a typical replicated series.
n_rep <- 4000
yrep <- lapply(c("B1", "B3", "B4"), function(m) {
  dr <- bayes$draws[[m]]
  dr <- dr[sample(nrow(dr), n_rep), ]
  mu <- mean_curve(m, dr, d$s)
  mu + matrix(rnorm(length(mu), 0, dr$sigma_1), nrow(mu))
})
names(yrep) <- c("B1", "B3", "B4")

# Test statistics T(y). Each targets a feature the models could miss:
recent <- annual$year >= 2008
ppc_stats <- list(
  "SD of summers" = sd,
  "Hottest summer" = max,
  "Lag-1 autocorrelation" = function(y) cor(y[-1], y[-length(y)]),
  "Mean 2008–2026 minus mean 1938–2007" = function(y) mean(y[recent]) - mean(y[!recent])
)

ppc_long <- bind_rows(lapply(names(yrep), function(m) {
  bind_rows(lapply(names(ppc_stats), function(stat) {
    f <- ppc_stats[[stat]]
    tibble(model = m, stat = stat, value = apply(yrep[[m]], 1, f), observed = f(d$y))
  }))
}))

# posterior predictive p-value: P(T(y_rep) >= T(y)). Values near 0 or 1 mean
# the model rarely produces data like ours for that feature.
ppc_pvalues <- ppc_long |>
  group_by(model, stat) |>
  summarise(observed = first(observed), p_value = mean(value >= observed), .groups = "drop")
print(ppc_pvalues)

p_ppc <- ggplot(ppc_long, aes(value)) +
  geom_histogram(bins = 40, fill = "grey75", colour = NA) +
  geom_vline(aes(xintercept = observed), colour = col_cp, linewidth = 0.9) +
  facet_grid(model ~ stat, scales = "free_x") +
  labs(x = "Test statistic in replicated data (orange line: observed)", y = "Count") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), axis.text.y = element_blank())
ggsave(file.path(fig_dir, "bayes_ppc_stats.png"), p_ppc, width = 10, height = 5.5, dpi = 150)

p_dens <- ppc_dens_overlay(d$y, yrep$B3[1:100, ]) +
  labs(x = "Mean summer TMAX (°C)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
ggsave(file.path(fig_dir, "bayes_ppc_density.png"), p_dens, width = 6, height = 3.5, dpi = 150)

# change-year posterior under each prior variant
p_sens <- bind_rows(lapply(names(sens_runs), function(v) {
  tibble(variant = v, gamma = sens_runs[[v]]$draws$gamma)
})) |>
  mutate(variant = factor(variant, levels = names(sens_runs))) |>
  count(variant, gamma) |>
  group_by(variant) |>
  mutate(prob = n / sum(n)) |>
  ggplot(aes(gamma, prob)) +
  geom_col(fill = col_cp, width = 0.8) +
  facet_wrap(~ variant, ncol = 1, strip.position = "top") +
  labs(x = "Last year before the change (gamma)", y = "Posterior probability") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "bayes_sensitivity_changeyear.png"), p_sens, width = 8, height = 7, dpi = 150)

# ---- 5. save results for the report ----------------------------

saveRDS(list(prior_predictive = prior_predictive, sensitivity = sensitivity,
             sens_convergence = sens_convergence, b1_ar_summary = b1_ar_summary,
             ppc_pvalues = ppc_pvalues),
        file.path(out_dir, "bayesian_sensitivity.rds"))
