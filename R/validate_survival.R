#' Bootstrap confidence interval for Harrell's C-index
#'
#' The concordance (C) index with a non-parametric percentile bootstrap
#' interval. `risk` is a per-subject risk score where **higher = worse
#' survival** (e.g. a Cox linear predictor); the interval is obtained by
#' resampling subjects with replacement and recomputing C.
#'
#' @param time Numeric follow-up time per subject.
#' @param event Event indicator (1/`TRUE` = event, 0/`FALSE` = censored).
#' @param risk Numeric risk score, higher = higher hazard / shorter survival.
#' @param n_boot Number of bootstrap resamples (default 1000).
#' @param conf Confidence level (default 0.95).
#' @param seed Random seed (default 1).
#'
#' @return A named numeric vector `c(cindex, lower, upper)`.
#' @export
#'
#' @examples
#' set.seed(1)
#' risk <- rnorm(200)
#' time <- rexp(200, rate = exp(0.8 * risk))
#' event <- rbinom(200, 1, 0.7)
#' cindex_ci(time, event, risk, n_boot = 200)
# Harrell's C for a risk score (higher = worse survival). Internal.
harrell_c <- function(time, event, risk) {
  survival::concordance(
    survival::Surv(time, event) ~ risk, reverse = TRUE
  )$concordance
}

cindex_ci <- function(time, event, risk, n_boot = 1000L, conf = 0.95, seed = 1L) {
  ok <- is.finite(time) & is.finite(as.numeric(event)) & is.finite(risk) & time > 0
  time <- time[ok]
  event <- as.integer(event[ok])
  risk <- risk[ok]
  if (length(time) < 3L) stop("need at least 3 subjects with complete data")

  cstat <- function(idx) harrell_c(time[idx], event[idx], risk[idx])
  est <- cstat(seq_along(time))

  set.seed(seed)
  n <- length(time)
  boot <- vapply(seq_len(n_boot), function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    tryCatch(cstat(idx), error = function(e) NA_real_)
  }, numeric(1L))

  a <- (1 - conf) / 2
  ci <- stats::quantile(boot, c(a, 1 - a), na.rm = TRUE)
  c(cindex = est, lower = unname(ci[1L]), upper = unname(ci[2L]))
}

#' Rigorously validate a survival signature
#'
#' Evaluates a candidate gene / feature signature the way Phase 4 of the roadmap
#' asks for: a bootstrapped C-index, **always compared against a clinical-only
#' baseline** (and the two combined), with optional train → external-cohort
#' validation and a calibration table at a chosen time horizon. Models are plain
#' Cox proportional-hazards fits (`survival::coxph`); the signature's linear
#' predictor is the risk score.
#'
#' @param data Training data frame: one row per subject, with `time`, `event`,
#'   the `predictors`, and any `clinical` columns.
#' @param predictors Character vector of column names forming the signature to
#'   validate (e.g. genes surfaced by [integrate_diablo()] / [integrate_mofa()]).
#' @param time,event Column names of the follow-up time and 0/1 event indicator.
#' @param clinical Optional character vector of clinical covariate columns used
#'   for the baseline model (e.g. age, stage). Strongly recommended — a signature
#'   that cannot beat clinical variables alone has not earned its complexity.
#' @param test Optional external/held-out data frame with the same columns; when
#'   supplied, models are fit on `data` and evaluated on `test` (the honest
#'   setting). When `NULL`, evaluation is in-sample and optimistic.
#' @param horizon Optional time point for a calibration table (same units as
#'   `time`); `NULL` skips calibration.
#' @param n_groups Number of risk groups for calibration (default 4).
#' @param n_boot,seed Passed to [cindex_ci()].
#'
#' @return An object of class `survival_validation`: a list with `cindex` (a data
#'   frame, one row per model — `signature`, `clinical`, `combined` — with
#'   `cindex`/`lower`/`upper`), `calibration` (a data frame or `NULL`), and
#'   `meta`. Has `print` and [plot_calibration()] methods.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 300
#' g1 <- rnorm(n); g2 <- rnorm(n); age <- rnorm(n, 60, 10)
#' lp <- 0.7 * g1 - 0.5 * g2 + 0.02 * age
#' d <- data.frame(
#'   time = rexp(n, rate = exp(lp - mean(lp))),
#'   event = rbinom(n, 1, 0.7), g1 = g1, g2 = g2, age = age
#' )
#' validate_survival(d, predictors = c("g1", "g2"), clinical = "age",
#'   horizon = 1, n_boot = 200)
#' }
validate_survival <- function(data, predictors, time = "time", event = "event",
                              clinical = NULL, test = NULL, horizon = NULL,
                              n_groups = 4L, n_boot = 1000L, seed = 1L) {
  if (!is.data.frame(data)) stop("`data` must be a data frame")
  if (!is.character(predictors) || length(predictors) == 0L) {
    stop("`predictors` must be a non-empty character vector of column names")
  }
  eval_df <- if (is.null(test)) data else test
  if (!is.data.frame(eval_df)) stop("`test` must be a data frame")

  need <- c(time, event, predictors, clinical)
  for (nm in c("data", "eval")) {
    df <- if (nm == "data") data else eval_df
    miss <- setdiff(need, colnames(df))
    if (length(miss) > 0L) {
      stop(nm, " is missing columns: ", paste(utils::head(miss, 5L), collapse = ", "))
    }
  }
  if (is.null(clinical)) {
    warning("no `clinical` baseline supplied; reporting the signature alone. ",
            "Phase 4 is about beating a clinical-only model — pass `clinical`.")
  }

  # Cox linear predictor: fit on `data`, score on `eval_df`.
  lp <- function(vars) {
    bt <- paste0("`", vars, "`")
    f <- stats::as.formula(sprintf(
      "survival::Surv(`%s`, `%s`) ~ %s", time, event, paste(bt, collapse = " + ")
    ))
    fit <- survival::coxph(f, data = data)
    as.numeric(stats::predict(fit, newdata = eval_df, type = "lp"))
  }

  et <- eval_df[[time]]
  ee <- eval_df[[event]]
  rows <- list(signature = cindex_ci(et, ee, lp(predictors), n_boot, seed = seed))
  if (!is.null(clinical)) {
    rows$clinical <- cindex_ci(et, ee, lp(clinical), n_boot, seed = seed)
    rows$combined <- cindex_ci(et, ee, lp(c(clinical, predictors)), n_boot, seed = seed)
  }
  cindex <- data.frame(
    model = names(rows),
    do.call(rbind, rows),
    row.names = NULL
  )

  calibration <- NULL
  if (!is.null(horizon)) {
    calibration <- calibrate_survival(
      data, eval_df, predictors, time, event, horizon, n_groups
    )
  }

  structure(
    list(
      cindex = cindex,
      calibration = calibration,
      meta = list(
        n_train = nrow(data), n_eval = nrow(eval_df),
        external = !is.null(test), horizon = horizon,
        n_predictors = length(predictors)
      )
    ),
    class = "survival_validation"
  )
}

# Calibration at a horizon: fit Cox on train, predict S(horizon) per eval
# subject, bin into risk groups, compare mean predicted vs observed KM survival.
calibrate_survival <- function(train, eval, predictors, time, event,
                               horizon, n_groups) {
  bt <- paste0("`", predictors, "`")
  f <- stats::as.formula(sprintf(
    "survival::Surv(`%s`, `%s`) ~ %s", time, event, paste(bt, collapse = " + ")
  ))
  fit <- survival::coxph(f, data = train)
  sf <- survival::survfit(fit, newdata = eval)
  pred <- as.numeric(summary(sf, times = horizon, extend = TRUE)$surv)

  qs <- stats::quantile(pred, probs = seq(0, 1, length.out = n_groups + 1L),
                        na.rm = TRUE)
  grp <- cut(pred, unique(qs), include.lowest = TRUE, labels = FALSE)

  parts <- lapply(sort(unique(grp[!is.na(grp)])), function(g) {
    idx <- which(grp == g)
    km <- survival::survfit(
      survival::Surv(eval[[time]][idx], eval[[event]][idx]) ~ 1
    )
    s <- summary(km, times = horizon, extend = TRUE)
    data.frame(
      group = g, n = length(idx),
      predicted = mean(pred[idx], na.rm = TRUE),
      observed = s$surv, obs_lower = s$lower, obs_upper = s$upper
    )
  })
  do.call(rbind, parts)
}

#' @export
print.survival_validation <- function(x, ...) {
  m <- x$meta
  cat(sprintf(
    "Survival validation | %s evaluation | train n=%d, eval n=%d | %d predictors\n",
    if (m$external) "EXTERNAL (held-out)" else "in-sample (optimistic)",
    m$n_train, m$n_eval, m$n_predictors
  ))
  ci <- x$cindex
  ci[, c("cindex", "lower", "upper")] <- round(ci[, c("cindex", "lower", "upper")], 3)
  cat("\nC-index (bootstrap 95% CI):\n")
  print(ci, row.names = FALSE)
  if (!is.null(x$calibration)) {
    cat(sprintf("\nCalibration at t=%s (predicted vs observed survival):\n", m$horizon))
    cal <- x$calibration
    cal[, c("predicted", "observed")] <- round(cal[, c("predicted", "observed")], 3)
    print(cal[, c("group", "n", "predicted", "observed")], row.names = FALSE)
  }
  invisible(x)
}

#' Plot a survival-validation calibration curve
#'
#' Predicted vs observed survival at the horizon, one point per risk group, with
#' the observed Kaplan–Meier confidence interval and the `y = x` line of perfect
#' calibration.
#'
#' @param validation A `survival_validation` object from [validate_survival()]
#'   built with a `horizon`.
#' @return A ggplot object.
#' @export
plot_calibration <- function(validation) {
  if (!inherits(validation, "survival_validation") ||
      is.null(validation$calibration)) {
    stop("`validation` must be a survival_validation built with a `horizon`")
  }
  cal <- validation$calibration
  ggplot2::ggplot(cal, ggplot2::aes(x = .data$predicted, y = .data$observed)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$obs_lower, ymax = .data$obs_upper),
      width = 0.02, colour = "grey40"
    ) +
    ggplot2::geom_point(size = 3, colour = "#2c7fb8") +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      x = sprintf("Predicted survival at t=%s", validation$meta$horizon),
      y = "Observed survival (Kaplan-Meier)",
      title = "Calibration"
    ) +
    ggplot2::theme_minimal()
}
