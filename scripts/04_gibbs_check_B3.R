# Section 8: our own Gibbs sampler for the jump model B3, as a check on mcp.
# Fits B3 with and without the slope. Run after scripts 01 and 02.

library(here)
library(coda)
source(here("R", "changepoint_tools.R"))

conv   <- readRDS(here("analysis", "results", "conventional_results.rds"))
annual <- conv$annual                   # year, n_days, y, s
y <- annual$y
s <- annual$s
n <- length(y)

# ---- settings -------------------------------------------------------------

m0 <- 35;  sd_level <- 3                # prior for int_1 and int_2
sd_slope <- 0.05                        # prior for s_1
sigma_sd <- 2                           # half-normal prior for sigma
k_years  <- 1953:2012                   # eligible first years of the new regime
split    <- 1987                        # boundary between the two modes (as in script 02)

n_iter   <- 10000                       # iterations per chain
burn     <- 1000
k_starts <- c(1955, 1970, 2000, 2010)   # one chain starts at each of these years
step_sd  <- 0.15                        # Metropolis step size for log(sigma)

# ---- the sampler ----------------------------------------------------------

gibbs_B3 <- function(slope, k_start) {
  # starting values
  k <- k_start
  int_1 <- mean(y); int_2 <- mean(y); s_1 <- 0; sigma <- sd(y)
  store <- matrix(0, nrow = n_iter, ncol = 5)
  colnames(store) <- c("int_1", "int_2", "s_1", "sigma", "k")

  for (it in 1:n_iter) {
    cp    <- k - 1982
    after <- annual$year >= k
    x     <- ifelse(after, s - cp, s)            # the slope's covariate

    # int_1 | rest: normal model on the years before the change
    z <- y[!after] - s_1 * x[!after]
    v <- 1 / (1 / sd_level^2 + length(z) / sigma^2)
    int_1 <- rnorm(1, v * (m0 / sd_level^2 + sum(z) / sigma^2), sqrt(v))

    # int_2 | rest: normal model on the years from the change on
    z <- y[after] - s_1 * x[after]
    v <- 1 / (1 / sd_level^2 + length(z) / sigma^2)
    int_2 <- rnorm(1, v * (m0 / sd_level^2 + sum(z) / sigma^2), sqrt(v))

    # s_1 | rest: regression through the origin on what the levels leave over
    if (slope) {
      r <- y - ifelse(after, int_2, int_1)
      v <- 1 / (1 / sd_slope^2 + sum(x^2) / sigma^2)
      s_1 <- rnorm(1, v * sum(x * r) / sigma^2, sqrt(v))
    }

    # sigma | rest: random-walk Metropolis on log(sigma), as in the Oct 9 code
    mu <- ifelse(after, int_2, int_1) + s_1 * x
    log_post <- function(sg) {
      sum(dnorm(y, mu, sg, log = TRUE)) + dnorm(sg, 0, sigma_sd, log = TRUE) + log(sg)
    }                                            # + log(sg): we step on the log scale
    proposed <- exp(log(sigma) + rnorm(1, 0, step_sd))
    if (log(runif(1)) < log_post(proposed) - log_post(sigma)) sigma <- proposed

    # k | rest: log-likelihood of every eligible change year, then sample one
    log_lik <- numeric(length(k_years))
    for (j in seq_along(k_years)) {
      aft <- annual$year >= k_years[j]
      mu_j <- ifelse(aft, int_2 + s_1 * (s - (k_years[j] - 1982)), int_1 + s_1 * s)
      log_lik[j] <- sum(dnorm(y, mu_j, sigma, log = TRUE))
    }
    k <- sample(k_years, 1, prob = exp(log_lik - max(log_lik)))

    store[it, ] <- c(int_1, int_2, s_1, sigma, k)
  }
  store[-(1:burn), ]
}

# ---- fit both versions ----------------------------------------------------

set.seed(638)
fits <- list()
for (version in c("mean-only", "slope")) {
  cat("\n==========", version, "==========\n")
  chains <- lapply(k_starts, function(k0) gibbs_B3(slope = (version == "slope"), k0))
  draws  <- do.call(rbind, chains)
  k      <- draws[, "k"]
  early  <- k < split

  # convergence (coda, as in the Chapter 9 class code)
  k_chains <- mcmc.list(lapply(chains, function(ch) mcmc(ch[, "k"])))
  cat("R-hat for k:", round(gelman.diag(k_chains)$psrf[1], 2),
      "  effective sample size for k:", round(effectiveSize(k_chains)), "\n")
  cat("P(change before", split, ") by chain:",
      round(sapply(chains, function(ch) mean(ch[, "k"] < split)), 2), "\n")
  cat("P(change before", split, ") pooled:", round(mean(early), 3), "\n")

  # jump and slope within each mode
  jump <- draws[, "int_2"] - (draws[, "int_1"] + draws[, "s_1"] * (k - 1982))
  for (mode in c("early", "late")) {
    keep <- if (mode == "early") early else !early
    if (sum(keep) < 100) next
    cat(sprintf("%s mode: jump %.2f C (%.2f, %.2f), slope %.2f C/decade\n", mode,
                mean(jump[keep]), quantile(jump[keep], 0.025), quantile(jump[keep], 0.975),
                10 * mean(draws[keep, "s_1"])))
  }
  fits[[version]] <- list(chains = chains, draws = draws)
}

# ---- compare with the exact grid posterior --------------------------------
# The team's cp_grid_posterior() integrates the coefficients out exactly.
# We evaluate it at the same whole years our sampler uses.

exact <- cp_grid_posterior("B3", annual,
                           prior_mean = c(int_1 = m0, int_2 = m0, s_1 = 0),
                           prior_sd   = c(int_1 = sd_level, int_2 = sd_level, s_1 = sd_slope),
                           sigma_sd   = sigma_sd, cp_grid = k_years - 1982)
cat("\nExact grid (slope version): P(change before", split, ") =",
    round(sum(exact$prob[exact$gamma + 1 < split]), 3), "\n")

# ---- figures --------------------------------------------------------------

fig_dir <- here("analysis", "results", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

plot_check <- function() {
  par(mfrow = c(2, 1), mar = c(4, 4, 2.5, 1))
  for (version in names(fits)) {
    k_prob <- table(factor(fits[[version]]$draws[, "k"], levels = k_years)) /
      nrow(fits[[version]]$draws)
    bp <- barplot(k_prob, names.arg = k_years - 1, col = "grey75", border = NA,
                  las = 2, cex.names = 0.6, ylab = "Posterior probability",
                  ylim = c(0, max(k_prob, exact$prob) * 1.15),
                  main = paste0(version, ": Gibbs sampler (bars)",
                                if (version == "slope") " vs exact grid (points)"))
    if (version == "slope") points(bp, exact$prob, pch = 19, col = "#eb6834", cex = 0.7)
  }
  par(mfrow = c(1, 1))
}
plot_traces <- function() {
  par(mfrow = c(2, 1), mar = c(4, 4, 2.5, 1))
  for (version in names(fits)) {
    ch <- fits[[version]]$chains
    plot(ch[[1]][, "k"] - 1, type = "l", ylim = range(k_years) - 1,
         xlab = "Iteration", ylab = "Change year (gamma)",
         main = paste0(version, ": change-year traces, 4 chains"))
    for (i in 2:4) lines(ch[[i]][, "k"] - 1, col = i)
    abline(h = split - 0.5, lty = 2)
  }
  par(mfrow = c(1, 1))
}

for (f in c("check", "traces")) {
  plot_fn <- get(paste0("plot_", f))
  if (interactive()) plot_fn()
  png(file.path(fig_dir, paste0("gibbs_", f, ".png")), width = 1000, height = 750, res = 120)
  plot_fn()
  dev.off()
}

saveRDS(list(fits = fits, exact = exact), here("analysis", "results", "gibbs_check.rds"))
