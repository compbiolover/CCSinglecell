# Null-rigor diagnostics -------------------------------------------------------
#
# When cv_validate_survival() reports that the genomic signature adds nothing
# beyond clinical (dC ~ 0), that number alone cannot tell you *why*: is there
# genuinely no signal, or is the cohort simply too small/censored to detect the
# signal that exists? These three functions turn a bare null into a defensible,
# quantified claim:
#
#   * assess_null()    - is the observed dC distinguishable from noise?
#                        (permutation p-value for the combined - clinical gain)
#   * power_curve()    - what dC could this cohort even have detected?
#                        (semi-synthetic signal injection -> detectable floor)
#   * learning_curve() - is the null power-limited or signal-limited?
#                        (out-of-fold C vs training-set size)
#
# All three drive the same honest penalized-CV engine as cv_validate_survival()
# so their verdicts are directly comparable to it.

# Stratified subsample: draw `frac` of the rows within each event class so the
# event rate is preserved. Returns row indices into a length-`n` vector.
.strat_subsample <- function(event, frac) {
  keep <- integer(0)
  for (cls in c(0L, 1L)) {
    i <- which(event == cls)
    if (length(i) == 0L) next
    take <- max(1L, round(length(i) * frac))
    keep <- c(keep, sample(i, take))
  }
  sort(keep)
}

#' Is a null genomic gain distinguishable from noise? (permutation test)
#'
#' [cv_validate_survival()] can report that adding a genomic signature to the
#' clinical baseline yields `dC ~ 0`. This function decides whether that gain is
#' any larger than what a *randomly aligned* genomic block of the same size and
#' correlation structure would produce — the difference between "no signal" and
#' "signal too small to matter here."
#'
#' The observed statistic is the out-of-fold `combined - clinical` C-index gain
#' from a single event-stratified penalized-CV pass (the same engine as
#' [cv_validate_survival()]). The null distribution is built by **permuting only
#' the rows of the genomic predictor block** `n_perm` times — this breaks the
#' genomics ↔ outcome link while leaving the clinical ↔ outcome signal (and the
#' fold assignment) intact, so the test isolates the *increment* genomics claims
#' to add. The one-sided empirical p-value is
#' `(1 + #{null >= observed}) / (1 + n_perm)`.
#'
#' @param data Data frame: one row per subject with `time`, `event`,
#'   `predictors`, and `clinical` columns.
#' @param predictors Character vector of genomic feature columns (the block
#'   whose incremental value is being tested).
#' @param clinical Character vector of clinical covariate columns (the baseline
#'   held fixed under the null). Required.
#' @param time,event Column names of follow-up time and 0/1 event indicator.
#' @param n_perm Number of predictor-block permutations (default 200).
#' @param n_folds,n_repeats,alpha,seed Passed to [cv_validate_survival()].
#'   `n_repeats` is used only for the reported full-strength interval; each
#'   permutation uses a single fixed-fold pass for speed and comparability.
#'
#' @return An object of class `null_assessment`: a list with `observed` (the
#'   single-pass gain), `observed_ci` (the full `n_repeats` gain with interval
#'   and win-rate, from [cv_validate_survival()]), `null` (the vector of
#'   permuted gains), `p_value`, and `meta`. Has a `print` method and
#'   [plot_null_assessment()].
#' @seealso [power_curve()], [learning_curve()], [cv_validate_survival()]
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 200
#' g <- matrix(rnorm(n * 15), n, dimnames = list(NULL, paste0("g", 1:15)))
#' age <- rnorm(n, 60, 10)
#' # genes are pure noise; only age drives survival
#' d <- data.frame(time = rexp(n, exp(0.04 * (age - 60))),
#'   event = rbinom(n, 1, 0.7), age = age, g)
#' assess_null(d, predictors = paste0("g", 1:15), clinical = "age",
#'   n_perm = 50, n_repeats = 5)
#' }
assess_null <- function(data, predictors, clinical, time = "time",
                        event = "event", n_perm = 200L, n_folds = 5L,
                        n_repeats = 10L, alpha = 0.5, seed = 1L) {
  if (!is.data.frame(data)) stop("`data` must be a data frame")
  if (missing(clinical) || length(clinical) == 0L) {
    stop("`clinical` is required — the permutation null is about the genomic ",
         "*increment* over a fixed clinical baseline")
  }
  need <- c(time, event, predictors, clinical)
  miss <- setdiff(need, colnames(data))
  if (length(miss) > 0L) {
    stop("data is missing columns: ", paste(utils::head(miss, 5L), collapse = ", "))
  }

  # single fixed-fold CV pass -> combined - clinical gain (the test statistic)
  stat <- function(df) {
    cv_validate_survival(df, predictors, clinical, time, event,
      n_folds = n_folds, n_repeats = 1L, alpha = alpha, seed = seed)$delta$mean
  }

  # reported full-strength gain (interval + win-rate) for context
  observed_full <- cv_validate_survival(data, predictors, clinical, time, event,
    n_folds = n_folds, n_repeats = n_repeats, alpha = alpha, seed = seed)$delta

  t_obs <- stat(data)

  rows <- seq_len(nrow(data))
  null <- vapply(seq_len(n_perm), function(b) {
    set.seed(seed + b)                       # only the permutation varies
    perm <- sample(rows)
    df <- data
    df[, predictors] <- data[perm, predictors, drop = FALSE]
    tryCatch(stat(df), error = function(e) NA_real_)
  }, numeric(1L))

  n_ok <- sum(is.finite(null))
  p <- (1 + sum(null >= t_obs, na.rm = TRUE)) / (1 + n_ok)

  structure(list(
    observed = t_obs,
    observed_ci = observed_full,
    null = null,
    p_value = p,
    meta = list(n = nrow(data), n_perm = n_perm, n_perm_ok = n_ok,
      n_predictors = length(predictors), n_clinical = length(clinical),
      n_folds = n_folds, n_repeats = n_repeats, alpha = alpha)
  ), class = "null_assessment")
}

#' @export
print.null_assessment <- function(x, ...) {
  m <- x$meta
  cat(sprintf(
    "Permutation null for genomic increment | n=%d | %d genomic vs %d clinical | %d perms\n",
    m$n, m$n_predictors, m$n_clinical, m$n_perm_ok
  ))
  ci <- x$observed_ci
  cat(sprintf(
    "\nObserved gain (combined - clinical): dC = %+.3f  [full: %+.3f, %+.3f], wins %d%%\n",
    x$observed, ci$lower, ci$upper, round(100 * ci$win_rate)
  ))
  cat(sprintf(
    "Permuted-block gains:   mean %+.3f, 95%% <= %+.3f\n",
    mean(x$null, na.rm = TRUE),
    stats::quantile(x$null, 0.95, na.rm = TRUE)
  ))
  cat(sprintf("Permutation p-value: %.3f\n", x$p_value))
  if (x$p_value > 0.05) {
    cat("=> observed gain is within the permutation null: no evidence the genomic\n",
        "   block adds more than a randomly aligned block of the same size.\n", sep = "")
  } else {
    cat("=> observed gain exceeds the permutation null: the genomic block carries\n",
        "   real incremental signal (p <= 0.05).\n", sep = "")
  }
  invisible(x)
}

#' Detectable-effect (power) curve for an incremental genomic signature
#'
#' A `dC ~ 0` result is only interpretable once you know what `dC` the cohort
#' could have detected in the first place. This runs a **semi-synthetic power
#' analysis**: it preserves the real clinical linear predictor, sample size and
#' event rate, then injects a single synthetic prognostic feature of known
#' strength into the outcome *and* into the predictor block, and measures how
#' often the [cv_validate_survival()] harness recovers a positive gain. Sweeping
#' the injected effect size traces the harness's power curve and pins down the
#' **minimum detectable effect** at a target power.
#'
#' For each `effect_size` and simulation, a standardized signal `s` is drawn;
#' survival times are regenerated from the hazard
#' `exp(clinical_lp + effect_size * s)` (so the clinical signal is retained and
#' `s` is the only genomic signal present), `s` is appended to `predictors`, and
#' the harness is run. Detection is `dC` lower bound `> 0`. The reported
#' `mde` is the smallest injected effect whose detection rate reaches
#' `power_target`, with `mde_delta` its mean recovered `dC`.
#'
#' @param data,predictors,clinical,time,event As in [assess_null()]. `clinical`
#'   is required (its fitted linear predictor is the retained baseline).
#' @param effect_sizes Numeric vector of injected effect sizes on the
#'   hazard-log scale (default `c(0, 0.1, 0.2, 0.35, 0.5)`). `0` estimates the
#'   false-positive rate.
#' @param n_sim Simulations per effect size (default 20).
#' @param power_target Target detection rate for the minimum detectable effect
#'   (default 0.8).
#' @param n_folds,n_repeats,alpha,seed Passed to [cv_validate_survival()].
#'
#' @return An object of class `power_curve`: a list with `curve` (a data frame:
#'   `effect`, `mean_delta`, `detect_rate`), `mde`, `mde_delta`,
#'   `power_target`, and `meta`. Has a `print` method and [plot_power_curve()].
#' @seealso [assess_null()], [learning_curve()]
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 200
#' age <- rnorm(n, 60, 10)
#' g <- matrix(rnorm(n * 15), n, dimnames = list(NULL, paste0("g", 1:15)))
#' d <- data.frame(time = rexp(n, exp(0.04 * (age - 60))),
#'   event = rbinom(n, 1, 0.7), age = age, g)
#' power_curve(d, predictors = paste0("g", 1:15), clinical = "age",
#'   effect_sizes = c(0, 0.3, 0.6), n_sim = 5, n_repeats = 5)
#' }
power_curve <- function(data, predictors, clinical, time = "time",
                        event = "event",
                        effect_sizes = c(0, 0.1, 0.2, 0.35, 0.5),
                        n_sim = 20L, power_target = 0.8, n_folds = 5L,
                        n_repeats = 5L, alpha = 0.5, seed = 1L) {
  if (!is.data.frame(data)) stop("`data` must be a data frame")
  if (missing(clinical) || length(clinical) == 0L) {
    stop("`clinical` is required — it supplies the retained baseline signal")
  }
  need <- c(time, event, predictors, clinical)
  miss <- setdiff(need, colnames(data))
  if (length(miss) > 0L) {
    stop("data is missing columns: ", paste(utils::head(miss, 5L), collapse = ", "))
  }
  if ("planted_signal" %in% colnames(data)) {
    stop("`data` already has a 'planted_signal' column; rename it")
  }

  cl_f <- stats::as.formula(sprintf(
    "survival::Surv(`%s`, `%s`) ~ %s", time, event,
    paste0("`", clinical, "`", collapse = " + ")
  ))
  lp_clin <- as.numeric(stats::predict(survival::coxph(cl_f, data = data),
    type = "lp"))
  lp_clin <- lp_clin - mean(lp_clin)
  n <- nrow(data)
  preds_aug <- c(predictors, "planted_signal")

  one <- function(eff, s) {
    set.seed(seed * 1000L + as.integer(round(eff * 100)) + s)
    sig <- as.numeric(scale(stats::rnorm(n)))
    lp <- lp_clin + eff * sig
    df <- data
    df[[time]] <- stats::rexp(n, rate = exp(lp - mean(lp)))
    df[["planted_signal"]] <- sig
    v <- tryCatch(
      cv_validate_survival(df, preds_aug, clinical, time, event,
        n_folds = n_folds, n_repeats = n_repeats, alpha = alpha, seed = seed),
      error = function(e) NULL
    )
    if (is.null(v)) return(c(delta = NA_real_, detect = NA_real_))
    c(delta = v$delta$mean, detect = as.numeric(v$delta$lower > 0))
  }

  curve <- lapply(effect_sizes, function(eff) {
    sims <- vapply(seq_len(n_sim), function(s) one(eff, s), numeric(2L))
    data.frame(effect = eff,
      mean_delta = mean(sims["delta", ], na.rm = TRUE),
      detect_rate = mean(sims["detect", ], na.rm = TRUE))
  })
  curve <- do.call(rbind, curve)

  hit <- which(curve$detect_rate >= power_target)
  mde <- if (length(hit)) curve$effect[min(hit)] else NA_real_
  mde_delta <- if (!is.na(mde)) curve$mean_delta[curve$effect == mde][1L] else NA_real_

  structure(list(
    curve = curve, mde = mde, mde_delta = mde_delta,
    power_target = power_target,
    meta = list(n = n, n_sim = n_sim, n_predictors = length(predictors),
      n_clinical = length(clinical), n_folds = n_folds, n_repeats = n_repeats,
      alpha = alpha)
  ), class = "power_curve")
}

#' @export
print.power_curve <- function(x, ...) {
  m <- x$meta
  cat(sprintf(
    "Detectable-effect analysis | n=%d | %d sims/effect | %d genomic + 1 planted vs %d clinical\n",
    m$n, m$n_sim, m$n_predictors, m$n_clinical
  ))
  cv <- x$curve
  cv$mean_delta <- round(cv$mean_delta, 3)
  cv$detect_rate <- round(cv$detect_rate, 2)
  cat("\nInjected effect -> recovered gain & detection rate:\n")
  print(cv, row.names = FALSE)
  if (is.na(x$mde)) {
    cat(sprintf(
      "\n=> even the largest injected effect never reached %d%% power: this cohort\n   is essentially blind to an incremental genomic signal.\n",
      round(100 * x$power_target)
    ))
  } else {
    cat(sprintf(
      "\nMinimum detectable effect (%d%% power): effect >= %.2f  (mean recovered dC ~ %+.3f)\n",
      round(100 * x$power_target), x$mde, x$mde_delta
    ))
    cat("=> a true incremental gain below this floor would be missed here.\n")
  }
  invisible(x)
}

#' Learning curve: is a null power-limited or signal-limited?
#'
#' Traces out-of-fold C-index against training-set size for the `signature`,
#' `clinical`, and `combined` models. If the `combined` curve is still climbing
#' at the full sample size, the cohort is **power-limited** — more data could yet
#' reveal a genomic gain. If it has plateaued at (or below) the clinical curve,
#' the null is **signal-limited** — more of the same data will not help. This is
#' usually the single most convincing figure that a `dC ~ 0` result is real.
#'
#' Each `fraction` is subsampled `n_subsample` times (event-stratified), one
#' single-pass [cv_validate_survival()] is run per subsample, and the per-model
#' C-indices are averaged. A least-squares slope of the mean `combined` C-index
#' against training size summarizes the trend.
#'
#' @param data,predictors,clinical,time,event As in [assess_null()].
#' @param fractions Numeric training-size fractions to evaluate (default
#'   `c(0.4, 0.55, 0.7, 0.85, 1.0)`).
#' @param n_subsample Subsamples per fraction (default 10).
#' @param n_folds,alpha,seed Passed to [cv_validate_survival()] (each pass uses
#'   `n_repeats = 1`).
#'
#' @return An object of class `learning_curve`: a list with `curve` (mean
#'   C-index per model per fraction), `raw` (per-subsample rows),
#'   `combined_slope` (mean combined C-index per subject added), and `meta`. Has
#'   a `print` method and [plot_learning_curve()].
#' @seealso [assess_null()], [power_curve()]
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 300
#' age <- rnorm(n, 60, 10)
#' g <- matrix(rnorm(n * 15), n, dimnames = list(NULL, paste0("g", 1:15)))
#' d <- data.frame(time = rexp(n, exp(0.04 * (age - 60))),
#'   event = rbinom(n, 1, 0.7), age = age, g)
#' learning_curve(d, predictors = paste0("g", 1:15), clinical = "age",
#'   fractions = c(0.5, 0.75, 1.0), n_subsample = 3)
#' }
learning_curve <- function(data, predictors, clinical, time = "time",
                           event = "event",
                           fractions = c(0.4, 0.55, 0.7, 0.85, 1.0),
                           n_subsample = 10L, n_folds = 5L, alpha = 0.5,
                           seed = 1L) {
  if (!is.data.frame(data)) stop("`data` must be a data frame")
  if (missing(clinical) || length(clinical) == 0L) {
    stop("`clinical` is required — the curve compares combined vs clinical")
  }
  need <- c(time, event, predictors, clinical)
  miss <- setdiff(need, colnames(data))
  if (length(miss) > 0L) {
    stop("data is missing columns: ", paste(utils::head(miss, 5L), collapse = ", "))
  }
  ee <- as.integer(data[[event]])

  rows <- list()
  for (f in fractions) {
    for (s in seq_len(n_subsample)) {
      set.seed(seed * 100L + as.integer(round(f * 100)) + s)
      idx <- .strat_subsample(ee, f)
      if (length(idx) < 3L * n_folds) next
      v <- tryCatch(
        cv_validate_survival(data[idx, , drop = FALSE], predictors, clinical,
          time, event, n_folds = n_folds, n_repeats = 1L, alpha = alpha,
          seed = seed),
        error = function(e) NULL
      )
      if (is.null(v)) next
      ci <- stats::setNames(v$cindex$cindex, v$cindex$model)
      rows[[length(rows) + 1L]] <- data.frame(
        fraction = f, n = length(idx),
        signature = unname(ci["signature"]),
        clinical = unname(ci["clinical"]),
        combined = unname(ci["combined"])
      )
    }
  }
  if (length(rows) == 0L) stop("no subsample was large enough for CV; lower fractions or n_folds")
  raw <- do.call(rbind, rows)

  agg <- stats::aggregate(cbind(n, signature, clinical, combined) ~ fraction,
    data = raw, FUN = mean)
  agg <- agg[order(agg$fraction), ]

  combined_slope <- if (nrow(agg) >= 2L) {
    unname(stats::coef(stats::lm(combined ~ n, data = agg))[2L])
  } else NA_real_

  structure(list(
    curve = agg, raw = raw, combined_slope = combined_slope,
    meta = list(n_full = nrow(data), n_subsample = n_subsample,
      n_predictors = length(predictors), n_clinical = length(clinical),
      n_folds = n_folds, alpha = alpha)
  ), class = "learning_curve")
}

#' @export
print.learning_curve <- function(x, ...) {
  m <- x$meta
  cat(sprintf(
    "Learning curve | up to n=%d | %d subsamples/point | %d genomic vs %d clinical\n",
    m$n_full, m$n_subsample, m$n_predictors, m$n_clinical
  ))
  cv <- x$curve
  cv[, c("signature", "clinical", "combined")] <-
    round(cv[, c("signature", "clinical", "combined")], 3)
  cv$n <- round(cv$n)
  cat("\nMean out-of-fold C-index vs training size:\n")
  print(cv[, c("fraction", "n", "signature", "clinical", "combined")],
    row.names = FALSE)
  cat(sprintf("\nCombined-C slope: %+.5f per subject added\n", x$combined_slope))
  if (is.finite(x$combined_slope) && x$combined_slope > 1e-4) {
    cat("=> combined discrimination is still rising with n: POWER-limited —\n",
        "   a larger cohort could yet expose a genomic gain.\n", sep = "")
  } else {
    cat("=> combined discrimination has plateaued: SIGNAL-limited —\n",
        "   more of the same data will not reveal a genomic gain.\n", sep = "")
  }
  invisible(x)
}
