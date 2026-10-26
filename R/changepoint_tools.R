# Helpers for the Bayesian change-point analysis (scripts 02 and 03).
#
# Models (s = year - 1982). The B-models are the Bayesian twins of the
# conventional M-models in scripts/01_conventional_analysis.R:
#   B0  mu = int_1                                   (no change)
#   B1  mu = int_1 + s_1*s                           (linear trend)
#   B2  mu = int_1 + s_1*s + s_1_E2*s^2              (quadratic trend)
#   B3  before cp: int_1 + s_1*s,  after: int_2 + s_1*(s - cp_1)   (jump, same slope)
#   B4  int_1 + s_1*min(s, cp_1) + s_2*(max(s, cp_1) - cp_1)       (slope change, joined)
#   B5  before cp: int_1 + s_1*s,  after: int_2 + s_2*(s - cp_1)   (jump + slope change)
# Parameter names are the ones mcp uses, so the same names appear in the
# priors, the posterior draws and these formulas.

#' mcp model formulas (lists of segments) for the six candidate models
#'
#' @param ar Add AR(1) errors (mcp's `ar(1)` term) to the first segment.
mcp_model <- function(model, ar = FALSE) {
  # mcp needs the intercept written out ("1 +"), so the formulas are spelled out
  first <- if (ar) {
    switch(model,
      B0 = y ~ 1 + ar(1),
      B1 = y ~ 1 + s + ar(1),
      B2 = y ~ 1 + s + I(s^2) + ar(1),
      y ~ 1 + s + ar(1))
  } else {
    switch(model,
      B0 = y ~ 1,
      B1 = y ~ 1 + s,
      B2 = y ~ 1 + s + I(s^2),
      y ~ 1 + s)
  }
  switch(model,
    B0 = , B1 = , B2 = list(first),
    B3 = list(first, ~ 1 + s),   # new intercept: the line can jump at cp_1
    B4 = list(first, ~ 0 + s),   # no new intercept: the line stays joined at cp_1
    B5 = list(first, ~ 1 + s))
}

#' Names of the regression coefficients in each model
coef_names <- function(model) {
  switch(model,
    B0 = "int_1",
    B1 = c("int_1", "s_1"),
    B2 = c("int_1", "s_1", "s_1_E2"),
    B3 = c("int_1", "int_2", "s_1"),
    B4 = c("int_1", "s_1", "s_2"),
    B5 = c("int_1", "int_2", "s_1", "s_2"))
}

#' Build the list of priors that mcp expects
#'
#' Every regression coefficient gets a normal prior N(prior_mean, prior_sd^2),
#' sigma_1 a half-normal prior with SD sigma_sd, and cp_1 a uniform prior on
#' (cp_lower, cp_upper). mcp writes priors in JAGS syntax, except that the
#' second argument of dnorm() is a standard deviation, not a precision.
#'
#' @param prior_mean,prior_sd Named numeric vectors (names as in coef_names()).
mcp_priors <- function(model, prior_mean, prior_sd, sigma_sd,
                       cp_lower = NULL, cp_upper = NULL, ar = FALSE) {
  priors <- list()
  for (p in coef_names(model)) {
    priors[[p]] <- sprintf("dnorm(%g, %g)", prior_mean[[p]], prior_sd[[p]])
  }
  priors$sigma_1 <- sprintf("dnorm(0, %g) T(0, )", sigma_sd)
  if (model %in% c("B3", "B4", "B5")) {
    priors$cp_1 <- sprintf("dunif(%g, %g)", cp_lower, cp_upper)
  }
  # B3: tie the slope after the change to the slope before it
  if (model == "B3") priors$s_2 <- "s_1"
  if (ar) priors$ar1_1 <- "dunif(-1, 1)"   # any stationary AR(1) coefficient equally likely
  priors
}

#' Reproducible JAGS initial values: one random-number seed per chain
chain_inits <- function(n_chains, seed) {
  lapply(seq_len(n_chains), function(i) {
    list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = seed + i)
  })
}

#' Fit one model with mcp
#'
#' @param mcmc List with chains, iter, adapt and seed.
#' @param extra_inits Optional list (one element per chain) of starting values.
fit_mcp <- function(model, data, priors, mcmc, ar = FALSE, extra_inits = NULL) {
  inits <- chain_inits(mcmc$chains, mcmc$seed)
  if (!is.null(extra_inits)) inits <- Map(c, inits, extra_inits)
  mcp::mcp(mcp_model(model, ar), data = data, prior = priors, par_x = "s",
           sample = "post", chains = mcmc$chains, iter = mcmc$iter,
           adapt = mcmc$adapt, cores = 1, inits = inits)
}

#' Posterior draws of an mcp fit as a data frame (one row per draw)
get_draws <- function(fit) as.data.frame(as.matrix(fit$mcmc_post))

#' Mean temperature mu(s) for every posterior draw
#'
#' @param draws Data frame of posterior draws (columns named as in mcp).
#' @param s Vector of time points (year - 1982).
#' @return A matrix with one row per draw and one column per time point.
mean_curve <- function(model, draws, s) {
  S <- matrix(s, nrow = nrow(draws), ncol = length(s), byrow = TRUE)
  if (model %in% c("B0", "B1", "B2")) {
    mu <- matrix(draws$int_1, nrow(draws), length(s))
    if (model != "B0") mu <- mu + draws$s_1 * S
    if (model == "B2") mu <- mu + draws$s_1_E2 * S^2
    return(mu)
  }
  cp <- draws$cp_1          # one value per row (draw), recycled across columns
  after <- S >= cp          # TRUE once the change point has passed
  switch(model,
    B3 = ifelse(after, draws$int_2 + draws$s_1 * (S - cp), draws$int_1 + draws$s_1 * S),
    B4 = draws$int_1 + draws$s_1 * pmin(S, cp) + draws$s_2 * (pmax(S, cp) - cp),
    B5 = ifelse(after, draws$int_2 + draws$s_2 * (S - cp), draws$int_1 + draws$s_1 * S))
}

#' Add the scientifically interesting quantities to every posterior draw
#'
#'   cp_year       change point on the calendar-year scale
#'   gamma         last year BEFORE the change (same definition as script 01,
#'                 where step = I(year > gamma)); the new segment starts at the
#'                 first s >= cp_1
#'   jump          size of the jump in level at the change point (deg C)
#'   trend_decade  warming rate in the first segment (deg C per decade)
#'   slope_change_decade  change in warming rate at the change point
#'   total_warming mean temperature in the last year minus the first year
#'
#' @param s_range First and last s in the data (c(-44, 44) for 1938-2026).
add_derived <- function(draws, model, s_range, mid_year = 1982) {
  if ("cp_1" %in% names(draws)) {
    draws$cp_year <- mid_year + draws$cp_1
    draws$gamma   <- mid_year + ceiling(draws$cp_1) - 1
  }
  if (model %in% c("B3", "B5")) {
    # segment 1 ends at int_1 + s_1*cp_1 and segment 2 starts at int_2
    draws$jump <- draws$int_2 - (draws$int_1 + draws$s_1 * draws$cp_1)
  }
  if (model %in% c("B4", "B5")) draws$slope_change_decade <- 10 * (draws$s_2 - draws$s_1)
  draws$trend_decade <- if ("s_1" %in% names(draws)) 10 * draws$s_1 else 0
  mu_ends <- mean_curve(model, draws, s_range)
  draws$total_warming <- mu_ends[, 2] - mu_ends[, 1]
  draws
}

#' Posterior mean, median and 95% equal-tailed credible interval
summ <- function(x) {
  tibble::tibble(mean = mean(x), median = stats::median(x),
                 lower = unname(stats::quantile(x, 0.025)),
                 upper = unname(stats::quantile(x, 0.975)))
}

#' Pointwise log-likelihood matrix (draws x years) for PSIS-LOO, iid errors
log_lik_matrix <- function(model, draws, data) {
  mu <- mean_curve(model, draws, data$s)
  Y  <- matrix(data$y, nrow(mu), ncol(mu), byrow = TRUE)
  stats::dnorm(Y, mu, draws$sigma_1, log = TRUE)
}

# ---- exact grid posterior of the change point ---------------------------
#
# For a FIXED change point cp (and fixed sigma, rho) every model above is an
# ordinary linear regression y = X beta + e with normal priors on beta. Then
# beta can be integrated out exactly:
#     y | cp, sigma ~ N(X m, sigma^2 I + X V X'),   m, V = prior mean, covariance.
# Evaluating this on a grid of cp and sigma values gives p(cp | y) without
# MCMC. We use it (1) to check the MCMC and (2) to weight MCMC runs that are
# stuck in different modes of a multimodal posterior.

#' Design matrix X for a change-point model at a fixed change point
cp_design <- function(model, s, cp) {
  after <- as.numeric(s >= cp)
  switch(model,
    B3 = cbind(int_1 = 1 - after, int_2 = after,
               s_1 = ifelse(after == 1, s - cp, s)),
    B4 = cbind(int_1 = 1, s_1 = pmin(s, cp), s_2 = pmax(s, cp) - cp),
    B5 = cbind(int_1 = 1 - after, int_2 = after,
               s_1 = (1 - after) * s, s_2 = after * (s - cp)))
}

#' log p(y | cp, sigma, rho), with the regression coefficients integrated out
#'
#' With AR(1) errors we follow mcp's (conditional) likelihood, which is the
#' same as quasi-differencing: y_t - rho*y_{t-1} = (X_t - rho*X_{t-1}) beta + e_t
#' for t >= 2, with the first observation left as it is.
#' The normal density is evaluated with the Woodbury identity, so only a small
#' (number of coefficients) x (number of coefficients) matrix is factorized.
cp_log_marginal <- function(y, X, prior_mean, prior_sd, sigma, rho = 0) {
  if (rho != 0) {
    n <- length(y)
    y <- c(y[1], y[-1] - rho * y[-n])
    X <- rbind(X[1, ], X[-1, , drop = FALSE] - rho * X[-n, , drop = FALSE])
  }
  n <- length(y)
  r <- y - X %*% prior_mean
  A <- diag(1 / prior_sd^2, ncol(X)) + crossprod(X) / sigma^2
  cA <- chol(A)
  b <- backsolve(cA, crossprod(X, r) / sigma^2, transpose = TRUE)
  quad   <- sum(r^2) / sigma^2 - sum(b^2)
  logdet <- 2 * n * log(sigma) + sum(log(prior_sd^2)) + 2 * sum(log(diag(cA)))
  -0.5 * (n * log(2 * pi) + logdet + quad)
}

log_sum_exp <- function(x) max(x) + log(sum(exp(x - max(x))))

#' Posterior of the change point on a grid
#'
#' The prior on cp is uniform over cp_grid, the prior on sigma is half-normal
#' with SD sigma_sd, and the prior on rho (if rho_grid has more than one value)
#' is uniform. sigma and rho are integrated out numerically on their grids.
#'
#' @return A tibble with columns cp, gamma (last year before the change, as in
#'   script 01) and prob (posterior probability of each grid point).
cp_grid_posterior <- function(model, data, prior_mean, prior_sd, sigma_sd,
                              cp_grid, sigma_grid = seq(0.7, 2, by = 0.02),
                              rho_grid = 0, mid_year = 1982) {
  m  <- prior_mean[coef_names(model)]
  sd <- prior_sd[coef_names(model)]
  sr <- expand.grid(sigma = sigma_grid, rho = rho_grid)
  log_post <- vapply(cp_grid, function(cp) {
    X <- cp_design(model, data$s, cp)
    lp <- mapply(function(sg, rho) cp_log_marginal(data$y, X, m, sd, sg, rho),
                 sr$sigma, sr$rho) +
      stats::dnorm(sr$sigma, 0, sigma_sd, log = TRUE)
    log_sum_exp(lp)
  }, numeric(1))
  prob <- exp(log_post - max(log_post))
  tibble::tibble(cp = cp_grid, gamma = mid_year + ceiling(cp_grid) - 1,
                 prob = prob / sum(prob))
}

#' Fit a change-point model whose posterior may have several modes
#'
#' JAGS updates one parameter at a time, so a chain that starts near one mode
#' of the change point rarely crosses to another mode. We therefore:
#'   1. compute the exact grid posterior of cp, and from it the posterior
#'      probability that cp is below / above `split`;
#'   2. fit the model with mcp twice, with the cp prior restricted to each side
#'      of `split` (each run then only has to explore one mode);
#'   3. combine the two sets of draws in proportion to the probabilities from 1.
#' This gives draws from the full posterior because
#'   p(theta | y) = P(early | y) p(theta | y, early) + P(late | y) p(theta | y, late).
#'
#' @return A list with the combined `draws` (column `region` says which run
#'   each draw came from), the two mcp `fits`, the `grid` posterior and the
#'   region `weights`.
fit_cp_by_region <- function(model, data, prior_mean, prior_sd, sigma_sd,
                             cp_lower, cp_upper, split, mcmc, ar = FALSE,
                             cp_step = 0.25, rho_grid = 0) {
  grid <- cp_grid_posterior(model, data, prior_mean, prior_sd, sigma_sd,
                            cp_grid = seq(cp_lower + cp_step / 2, cp_upper, by = cp_step),
                            rho_grid = if (ar) rho_grid else 0)
  weights <- c(early = sum(grid$prob[grid$cp < split]),
               late  = sum(grid$prob[grid$cp >= split]))

  bounds <- list(early = c(cp_lower, split), late = c(split, cp_upper))
  fits <- lapply(names(bounds), function(region) {
    priors <- mcp_priors(model, prior_mean, prior_sd, sigma_sd,
                         bounds[[region]][1], bounds[[region]][2], ar = ar)
    fit_mcp(model, data, priors, mcmc, ar = ar)
  })
  names(fits) <- names(bounds)

  # combine: take round(weight * total) draws from each run
  n_total <- mcmc$chains * mcmc$iter
  n_early <- round(weights[["early"]] * n_total)
  take <- c(early = n_early, late = n_total - n_early)
  draws <- dplyr::bind_rows(lapply(names(fits), function(region) {
    d <- get_draws(fits[[region]])
    d[sample(nrow(d), take[[region]]), , drop = FALSE] |>
      dplyr::mutate(region = region)
  }))

  list(draws = draws, fits = fits, grid = grid, weights = weights)
}
