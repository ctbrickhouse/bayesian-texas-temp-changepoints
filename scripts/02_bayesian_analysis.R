library(here)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(mcp)
library(loo)
library(posterior)

source(here("R", "changepoint_tools.R"))

# Bayesian version of the conventional analysis in 01_conventional_analysis.R.
# Run 01_conventional_analysis.R first: this script re-uses its Y(t) series
# and its results so the two analyses are directly comparable.
#
# Candidate mean models (s = year - 1982, the middle year). Each B-model is the
# Bayesian twin of the conventional M-model with the same number:
#   B0  no change          B3  jump in level (same slope before and after)
#   B1  linear trend       B4  change in slope (lines joined at the change)
#   B2  quadratic trend    B5  jump + change in slope
# The exact equations are at the top of R/changepoint_tools.R.
# All models use iid normal errors, y ~ N(mu, sigma^2), because the
# conventional analysis found no residual autocorrelation once a change point
# is in the model (AR(1) errors are checked in 03_bayesian_sensitivity.R).
#
# The models are fit with the mcp package, which writes the JAGS code for us
# and samples the posterior with MCMC (Gibbs sampling).

set.seed(638)

# MCMC settings (used for every mcp fit)
mcmc <- list(
  chains = 4,
  iter   = 20000,   # kept draws per chain
  adapt  = 3000,    # adaptation + burn-in draws per chain (discarded)
  seed   = 638      # JAGS random-number seeds are 639, 640, ... (one per chain)
)

out_dir <- here("analysis", "results")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

col_data  <- "grey55"
col_m1    <- "#2a78d6"
col_cp    <- "#eb6834"
col_early <- "#8a5cc2"
col_b4    <- "#2e9e6b"

# ---- 0. load Y(t) and the conventional results -----------------

conv <- readRDS(file.path(out_dir, "conventional_results.rds"))
annual <- conv$annual            # year, n_days, y, s
mid_year <- 1982                 # s = year - mid_year, as in script 01

# mcp wants a data frame that contains only the variables in the formulas
d <- annual |> select(s, y) |> as.data.frame()

# ---- 1. priors -------------------------------------------------

# Normal priors N(mean, sd^2) for the regression coefficients:
#   int_1, int_2  mean summer TMAX (deg C) in 1982 / just after the change:
#                 Austin summers are about 35 deg C; SD 3 allows anything
#                 from about 29 to 41 deg C
#   s_1, s_2      slope in deg C per year: SD 0.05 puts 95% of the prior
#                 within +/- 1 deg C per decade, about five times the
#                 observed rate of global warming
#   s_1_E2        curvature of the quadratic trend
prior_mean <- c(int_1 = 35, int_2 = 35, s_1 = 0,    s_2 = 0,    s_1_E2 = 0)
prior_sd   <- c(int_1 = 3,  int_2 = 3,  s_1 = 0.05, s_2 = 0.05, s_1_E2 = 0.002)
sigma_sd   <- 2      # half-normal prior on the year-to-year SD sigma_1
cp_lower   <- -30    # change point uniform on 1952-2012: the same 15%
cp_upper   <- 30     #   trimming of the record as gamma_grid in script 01

# ---- 2. models without a change point: B0, B1, B2 --------------

fits <- list()
for (m in c("B0", "B1", "B2")) {
  cat("\n==== fitting", m, "====\n")
  fits[[m]] <- fit_mcp(m, d, mcp_priors(m, prior_mean, prior_sd, sigma_sd), mcmc)
}
draws <- lapply(fits, get_draws)

# ---- 3. change-point models: B3, B4, B5 ------------------------

# 3a. A first, naive fit of B3 shows the problem. We start the four chains at
#     change points spread over the whole prior range (1957, 1974, 1990, 2007).
naive_inits <- lapply(c(-25, -8, 8, 25), function(cp) list(cp_1 = cp))
naive_fit <- fit_mcp("B3", d,
                     mcp_priors("B3", prior_mean, prior_sd, sigma_sd, cp_lower, cp_upper),
                     mcmc, extra_inits = naive_inits)
naive_check <- summarise_draws(as_draws_df(naive_fit$mcmc_post), "mean", "rhat", "ess_bulk")
naive_cp_by_chain <- sapply(naive_fit$mcmc_post, function(ch) mean(ch[, "cp_1"]) + mid_year)
print(naive_check)
cat("Mean change year in each chain:", round(naive_cp_by_chain, 1), "\n")
# Chains that start in the 1960s stay near the 1960s and chains that start in
# the 2000s stay near 2007, so R-hat is far above 1.01: the posterior of the
# change point has two separate modes and the Gibbs sampler cannot move
# between them.

# 3b. Fix (see fit_cp_by_region() in R/changepoint_tools.R): compute the exact
#     posterior probability of each mode on a grid, run mcp separately on each
#     side of the gap between the modes, and combine the two runs with those
#     probabilities. The split at s = 5 (1987) lies between the two modes.
split <- 5
cp_fits <- list()
for (m in c("B3", "B4", "B5")) {
  cat("\n==== fitting", m, "(two regions) ====\n")
  cp_fits[[m]] <- fit_cp_by_region(m, d, prior_mean, prior_sd, sigma_sd,
                                   cp_lower, cp_upper, split, mcmc)
  draws[[m]] <- cp_fits[[m]]$draws
  print(cp_fits[[m]]$weights)
}

# ---- 4. convergence diagnostics --------------------------------

# For every parameter of every MCMC run: rank-normalized R-hat (should be
# < 1.01) and bulk/tail effective sample sizes (should be well above 400).
all_runs <- c(fits,
              setNames(lapply(cp_fits, `[[`, "fits") |> unlist(recursive = FALSE),
                       paste(rep(names(cp_fits), each = 2), c("early", "late"))))
convergence <- bind_rows(lapply(names(all_runs), function(run) {
  summarise_draws(as_draws_df(all_runs[[run]]$mcmc_post),
                  "mean", "sd", "rhat", "ess_bulk", "ess_tail") |>
    filter(variable != "s_2" | !startsWith(run, "B3")) |>   # in B3, s_2 is a copy of s_1
    mutate(run = run, .before = 1)
}))
print(convergence, n = Inf)
cat("Largest R-hat:", round(max(convergence$rhat), 4),
    " smallest bulk ESS:", round(min(convergence$ess_bulk)),
    " smallest tail ESS:", round(min(convergence$ess_tail)), "\n")

# Check the MCMC against the exact grid posterior: within each region, the
# MCMC distribution of the change year should match the grid.
grid_quantile <- function(x, p, q) x[which(cumsum(p) / sum(p) >= q)[1]]
mcmc_vs_grid <- bind_rows(lapply(names(cp_fits), function(m) {
  bind_rows(lapply(c("early", "late"), function(region) {
    g <- cp_fits[[m]]$grid |> filter(if (region == "early") cp < split else cp >= split)
    cp_mcmc <- get_draws(cp_fits[[m]]$fits[[region]])$cp_1 + mid_year
    tibble(model = m, region = region,
           weight = cp_fits[[m]]$weights[[region]],
           mcmc_mean = mean(cp_mcmc),
           grid_mean = sum((g$cp + mid_year) * g$prob) / sum(g$prob),
           mcmc_lower = unname(quantile(cp_mcmc, 0.025)),
           grid_lower = grid_quantile(g$cp + mid_year, g$prob, 0.025),
           mcmc_upper = unname(quantile(cp_mcmc, 0.975)),
           grid_upper = grid_quantile(g$cp + mid_year, g$prob, 0.975))
  }))
}))
print(mcmc_vs_grid, width = Inf)

# ---- 5. model comparison with leave-one-out cross-validation ---

# elpd_loo = expected log predictive density for a new summer, estimated with
# Pareto-smoothed importance sampling (PSIS-LOO). Higher is better.
# r_eff = 1 treats the draws as independent; it only affects the Monte Carlo
# error of the estimate, not the estimate itself.
loos <- lapply(names(draws), function(m) {
  loo(log_lik_matrix(m, draws[[m]], d), r_eff = rep(1, nrow(d)))
})
names(loos) <- names(draws)
loo_table <- loo_compare(loos)
print(loo_table)

# turn the comparison into model weights (how much each model contributes to
# the best predictive average of the six models)
loo_weights <- tibble(
  model = names(loos),
  elpd_loo = sapply(loos, function(l) l$estimates["elpd_loo", "Estimate"]),
  elpd_se  = sapply(loos, function(l) l$estimates["elpd_loo", "SE"]),
  p_loo    = sapply(loos, function(l) l$estimates["p_loo", "Estimate"]),
  max_pareto_k = sapply(loos, function(l) max(l$diagnostics$pareto_k)),
  stacking   = as.numeric(loo_model_weights(loos, method = "stacking")),
  pseudo_bma = as.numeric(loo_model_weights(loos, method = "pseudobma"))
) |>
  left_join(as.data.frame(loo_table) |> select(model, elpd_diff, se_diff), by = "model") |>
  arrange(desc(elpd_loo))
print(loo_weights, width = Inf)

# ---- 6. posterior summaries ------------------------------------

# Add the scientifically interesting quantities (change year, jump size,
# warming rate, total warming 1938-2026, ...) to every draw; see add_derived()
# in R/changepoint_tools.R for the definitions.
s_range <- range(d$s)
for (m in names(draws)) draws[[m]] <- add_derived(draws[[m]], m, s_range, mid_year)

key_quantities <- c("cp_year", "jump", "trend_decade", "slope_change_decade",
                    "total_warming", "sigma_1")
summarise_model <- function(dr, model) {
  qs <- intersect(key_quantities, names(dr))
  bind_rows(lapply(qs, function(q) summ(dr[[q]]) |> mutate(quantity = q, .before = 1)))
}

# whole posterior of each model ...
posterior_summary <- bind_rows(lapply(names(draws), function(m) {
  summarise_model(draws[[m]], m) |> mutate(model = m, region = "all", .before = 1)
}))
# ... and, for the change-point models, separately within each mode
posterior_by_region <- bind_rows(lapply(c("B3", "B4", "B5"), function(m) {
  bind_rows(lapply(c("early", "late"), function(r) {
    summarise_model(filter(draws[[m]], region == r), m) |>
      mutate(model = m, region = r, weight = cp_fits[[m]]$weights[[r]], .before = 1)
  }))
}))
print(posterior_summary, n = Inf, width = Inf)
print(posterior_by_region, n = Inf, width = Inf)

best_cp_loo <- loo_weights |> filter(model %in% c("B3", "B4", "B5")) |> slice(1) |> pull(model)
cat("Best change-point model by LOO:", best_cp_loo, "\n")

# The rest of the script focuses on B3 (an abrupt jump in level), the Bayesian
# twin of M3, the model the conventional analysis selected. B4 (a gradual
# change in slope) is shown next to it because LOO cannot separate the two.
focus_cp <- "B3"
pb <- draws[[focus_cp]]

# posterior probabilities that answer the scientific question directly
post_probs <- tibble(
  statement = c(
    "Linear trend (B1): P(warming rate > 0)",
    "Linear trend (B1): P(warming rate > 0.2 °C/decade)",
    "Jump model (B3): P(mean 2026 > mean 1938)",
    "Jump model (B3): P(total warming > 1 °C)",
    "Jump model (B3): P(jump in 2000–2012)",
    "Jump model (B3): P(jump in 1952–1975)",
    "Jump model (B3): P(level jumps up, not down)",
    "Slope-change model (B4): P(warming rate increases)",
    "Slope-change model (B4): P(change in 1965–1990)"
  ),
  probability = c(
    mean(draws$B1$s_1 > 0),
    mean(draws$B1$trend_decade > 0.2),
    mean(pb$total_warming > 0),
    mean(pb$total_warming > 1),
    mean(pb$gamma >= 2000 & pb$gamma <= 2012),
    mean(pb$gamma >= 1952 & pb$gamma <= 1975),
    mean(pb$jump > 0),
    mean(draws$B4$slope_change_decade > 0),
    mean(draws$B4$gamma >= 1965 & draws$B4$gamma <= 1990)
  )
)
print(post_probs)

# posterior distribution of the change year (one row per year)
gamma_post <- pb |>
  count(gamma) |>
  mutate(prob = n / sum(n)) |>
  arrange(desc(prob))
print(head(gamma_post, 10))

# ---- 7. prediction for next summer (2027) ----------------------

# posterior predictive distribution of Y(2027): for every posterior draw,
# compute the mean for 2027 and add new year-to-year noise. This combines the
# uncertainty in all parameters (including the change year) with the noise.
predict_summer <- function(model, year) {
  mu <- mean_curve(model, draws[[model]], year - mid_year)[, 1]
  mu + rnorm(length(mu), 0, draws[[model]]$sigma_1)
}
pred_2027 <- bind_rows(lapply(c("B1", "B3", "B4"), function(m) {
  y_new <- predict_summer(m, 2027)
  summ(y_new) |> mutate(model = m, p_above_36 = mean(y_new > 36), .before = 1)
}))
print(pred_2027)

# ---- 8. figures ------------------------------------------------

# (a) data with posterior mean curve and 95% credible band: B1, B4 and the two
#     modes of B3
band <- function(model, dr, label) {
  mu <- mean_curve(model, dr[sample(nrow(dr), min(4000, nrow(dr))), ], d$s)
  tibble(year = annual$year, fit = label, mean = colMeans(mu),
         lower = apply(mu, 2, quantile, 0.025), upper = apply(mu, 2, quantile, 0.975))
}
w <- cp_fits[[focus_cp]]$weights
bands <- bind_rows(
  band("B1", draws$B1, "B1 linear trend"),
  band("B4", draws$B4, "B4 slope change"),
  band(focus_cp, filter(pb, region == "late"),
       sprintf("%s, late change (prob. %.2f)", focus_cp, w[["late"]])),
  band(focus_cp, filter(pb, region == "early"),
       sprintf("%s, early change (prob. %.2f)", focus_cp, w[["early"]]))
)
p_fit <- ggplot(bands, aes(year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = fit), alpha = 0.15) +
  geom_point(data = annual, aes(year, y), colour = col_data, size = 1.8) +
  geom_line(aes(y = mean, colour = fit), linewidth = 0.9) +
  scale_colour_manual(values = c(col_m1, col_early, col_cp, col_b4), name = NULL) +
  scale_fill_manual(values = c(col_m1, col_early, col_cp, col_b4), name = NULL) +
  labs(x = NULL, y = "Mean summer TMAX (°C)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))
ggsave(file.path(fig_dir, "bayes_fits.png"), p_fit, width = 8, height = 4.8, dpi = 150)

# (b) posterior of the change year next to the conventional bootstrap distribution
gamma_compare <- bind_rows(
  tibble(gamma = conv$boot_gamma, method = "1. Conventional parametric bootstrap, jump model (M3)"),
  tibble(gamma = pb$gamma, method = "2. Bayesian posterior, jump model (B3)"),
  tibble(gamma = draws$B4$gamma, method = "3. Bayesian posterior, slope-change model (B4)")
) |>
  count(method, gamma) |>
  group_by(method) |>
  mutate(prob = n / sum(n)) |>
  ungroup()
p_gamma <- ggplot(gamma_compare, aes(gamma, prob)) +
  geom_col(fill = col_cp, width = 0.8) +
  facet_wrap(~ method, ncol = 1) +
  labs(x = "Last year before the change (gamma)", y = "Probability") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "bayes_changeyear.png"), p_gamma, width = 8, height = 6, dpi = 150)

# (c) MCMC check: change-year probabilities from the combined MCMC draws (bars)
#     against the exact grid posterior (points), for B3, B4 and B5
check_df <- bind_rows(lapply(names(cp_fits), function(m) {
  mc <- draws[[m]] |> count(gamma) |> mutate(prob = n / sum(n), source = "MCMC (mcp)")
  gr <- cp_fits[[m]]$grid |> group_by(gamma) |> summarise(prob = sum(prob)) |>
    mutate(source = "Exact grid")
  bind_rows(mc, gr) |> mutate(model = m)
}))
p_check <- ggplot(check_df, aes(gamma, prob)) +
  geom_col(data = filter(check_df, source == "MCMC (mcp)"), fill = "grey75", width = 0.8) +
  geom_point(data = filter(check_df, source == "Exact grid"), colour = col_cp, size = 1.2) +
  facet_wrap(~ model, ncol = 1) +
  labs(x = "Last year before the change (gamma)", y = "Posterior probability") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "bayes_mcmc_vs_grid.png"), p_check, width = 8, height = 5.5, dpi = 150)

# (d) trace plots of the best change-point model, one panel per region
trace_pars <- intersect(c("cp_1", "int_1", "int_2", "s_1", "s_2", "sigma_1"),
                        colnames(as.matrix(cp_fits[[focus_cp]]$fits$late$mcmc_post)))
if (focus_cp == "B3") trace_pars <- setdiff(trace_pars, "s_2")
p_trace <- (plot_pars(cp_fits[[focus_cp]]$fits$early, pars = trace_pars, type = "trace", ncol = 2) +
              ggtitle("Early-change run")) /
  (plot_pars(cp_fits[[focus_cp]]$fits$late, pars = trace_pars, type = "trace", ncol = 2) +
     ggtitle("Late-change run"))
ggsave(file.path(fig_dir, "bayes_trace.png"), p_trace, width = 8, height = 9, dpi = 150)

# (e) prior vs posterior for the change year, the jump and the total warming
n_prior <- 40000
prior_draws <- tibble(
  cp_1  = runif(n_prior, cp_lower, cp_upper),
  int_1 = rnorm(n_prior, prior_mean[["int_1"]], prior_sd[["int_1"]]),
  int_2 = rnorm(n_prior, prior_mean[["int_2"]], prior_sd[["int_2"]]),
  s_1   = rnorm(n_prior, prior_mean[["s_1"]], prior_sd[["s_1"]]),
  s_2   = if (focus_cp == "B3") s_1 else rnorm(n_prior, prior_mean[["s_2"]], prior_sd[["s_2"]]),
  sigma_1 = abs(rnorm(n_prior, 0, sigma_sd))
) |> add_derived(focus_cp, s_range, mid_year)
pp_quantities <- c("Change year" = "cp_year", "Jump size (°C)" = "jump",
                   "Total warming 1938–2026 (°C)" = "total_warming")
pp_quantities <- pp_quantities[pp_quantities %in% names(pb)]
pp_long <- bind_rows(lapply(names(pp_quantities), function(lab) {
  bind_rows(tibble(quantity = lab, value = pb[[pp_quantities[[lab]]]], source = "Posterior"),
            tibble(quantity = lab, value = prior_draws[[pp_quantities[[lab]]]], source = "Prior"))
}))
p_pp <- ggplot(pp_long, aes(value, fill = source)) +
  geom_density(alpha = 0.5, colour = NA, adjust = 1.2) +
  facet_wrap(~ quantity, scales = "free") +
  scale_fill_manual(values = c(Posterior = col_cp, Prior = "grey70"), name = NULL) +
  labs(x = NULL, y = "Density") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "bayes_prior_posterior.png"), p_pp, width = 9, height = 3.5, dpi = 150)

# ---- 9. save results for the report ----------------------------

results <- list(
  settings = list(mcmc = mcmc, prior_mean = prior_mean, prior_sd = prior_sd,
                  sigma_sd = sigma_sd, cp_lower = cp_lower, cp_upper = cp_upper,
                  split = split, mid_year = mid_year),
  naive_check = naive_check, naive_cp_by_chain = naive_cp_by_chain,
  convergence = convergence, mcmc_vs_grid = mcmc_vs_grid,
  loo_weights = loo_weights,
  posterior_summary = posterior_summary, posterior_by_region = posterior_by_region,
  region_weights = lapply(cp_fits, `[[`, "weights"),
  post_probs = post_probs, gamma_post = gamma_post, pred_2027 = pred_2027,
  best_cp_loo = best_cp_loo, focus_cp = focus_cp, draws = draws
)
saveRDS(results, file.path(out_dir, "bayesian_results.rds"))
