#' Repeated cross-validated survival validation (penalized)
#'
#' The heavier-duty sibling of [validate_survival()]. A single train/test split
#' has high variance and unpenalized Cox overfits a wide signature; this runs
#' **repeated k-fold cross-validation** with **penalized (elastic-net) Cox** so
#' the signature is selected-and-shrunk honestly inside each training fold, and
#' answers the question that matters: *does adding the genomic signature to the
#' clinical variables improve out-of-fold discrimination, and how often?*
#'
#' Each repeat assigns event-stratified folds; every subject gets an
#' out-of-fold linear predictor from three models — the genomic `signature`
#' (elastic-net Cox over `predictors`), the `clinical` baseline (unpenalized
#' Cox), and `combined` (elastic-net Cox with the clinical terms forced in via a
#' zero penalty factor). One C-index per model per repeat gives a distribution,
#' not a single fragile number.
#'
#' @param data Data frame: one row per subject with `time`, `event`,
#'   `predictors`, and `clinical` columns.
#' @param predictors Character vector of candidate genomic feature columns
#'   (genes and/or miRNAs). Elastic-net selects among them within each fold.
#' @param clinical Character vector of clinical covariate columns (the baseline
#'   to beat). Required — the whole point is the comparison.
#' @param time,event Column names of follow-up time and 0/1 event indicator.
#' @param n_folds,n_repeats Cross-validation folds (default 5) and repeats
#'   (default 10).
#' @param alpha Elastic-net mixing parameter for the penalized fits (default
#'   0.5; 1 = lasso, 0 = ridge).
#' @param seed Random seed (default 1).
#'
#' @return An object of class `cv_survival_validation`: a list with `cindex` (a
#'   data frame — `signature`, `clinical`, `combined` — each with the mean C and
#'   a percentile interval across repeats), `delta` (the paired
#'   `combined - clinical` C-index difference: mean, interval, and the fraction
#'   of repeats where combined wins), `per_repeat` (the raw C matrix), and
#'   `meta`. Has a `print` method.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 250
#' g <- matrix(rnorm(n * 20), n, dimnames = list(NULL, paste0("g", 1:20)))
#' age <- rnorm(n, 60, 10)
#' lp <- 0.8 * g[, 1] - 0.6 * g[, 2] + 0.03 * age
#' d <- data.frame(time = rexp(n, exp(lp - mean(lp))), event = rbinom(n, 1, 0.7),
#'   age = age, g)
#' cv_validate_survival(d, predictors = paste0("g", 1:20), clinical = "age",
#'   n_repeats = 5)
#' }
cv_validate_survival <- function(data, predictors, clinical, time = "time",
                                 event = "event", n_folds = 5L, n_repeats = 10L,
                                 alpha = 0.5, seed = 1L) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("cv_validate_survival() requires the 'glmnet' package.")
  }
  if (!is.data.frame(data)) stop("`data` must be a data frame")
  if (missing(clinical) || length(clinical) == 0L) {
    stop("`clinical` is required — this function exists to test the signature ",
         "against a clinical baseline")
  }
  need <- c(time, event, predictors, clinical)
  miss <- setdiff(need, colnames(data))
  if (length(miss) > 0L) {
    stop("data is missing columns: ", paste(utils::head(miss, 5L), collapse = ", "))
  }

  data <- data[stats::complete.cases(data[, need, drop = FALSE]), , drop = FALSE]
  tt <- data[[time]]
  ee <- as.integer(data[[event]])
  keep <- is.finite(tt) & tt > 0 & !is.na(ee)
  data <- data[keep, , drop = FALSE]; tt <- tt[keep]; ee <- ee[keep]
  n <- nrow(data)
  if (n < 3L * n_folds) stop("too few complete cases for ", n_folds, "-fold CV")

  Xp <- as.matrix(data[, predictors, drop = FALSE])
  Xc <- as.matrix(data[, clinical, drop = FALSE])
  y <- survival::Surv(tt, ee)
  pen_combined <- c(rep(0, length(clinical)), rep(1, length(predictors)))

  # penalized Cox linear predictor for a fold (auto lambda via inner cv.glmnet)
  pen_lp <- function(x_tr, y_tr, x_te, penalty = NULL) {
    cv <- suppressMessages(glmnet::cv.glmnet(x_tr, y_tr, family = "cox",
      alpha = alpha, nfolds = 5L, cox.ties = "efron",
      penalty.factor = penalty %||% rep(1, ncol(x_tr))))
    as.numeric(stats::predict(cv, newx = x_te, s = "lambda.min"))
  }
  make_folds <- function(rng_seed) {
    set.seed(rng_seed)
    fold <- integer(n)
    for (cls in c(0L, 1L)) {                     # stratify by event
      i <- which(ee == cls)
      fold[i] <- sample(rep(seq_len(n_folds), length.out = length(i)))
    }
    fold
  }

  Cmat <- matrix(NA_real_, n_repeats, 3L,
    dimnames = list(NULL, c("signature", "clinical", "combined")))
  for (r in seq_len(n_repeats)) {
    fold <- make_folds(seed + r)
    lp_sig <- lp_clin <- lp_comb <- rep(NA_real_, n)
    for (k in seq_len(n_folds)) {
      te <- which(fold == k); tr_ <- which(fold != k)
      y_tr <- y[tr_]
      lp_sig[te]  <- pen_lp(Xp[tr_, , drop = FALSE], y_tr, Xp[te, , drop = FALSE])
      lp_comb[te] <- pen_lp(cbind(Xc, Xp)[tr_, , drop = FALSE], y_tr,
                            cbind(Xc, Xp)[te, , drop = FALSE], pen_combined)
      cfit <- survival::coxph(y_tr ~ Xc[tr_, , drop = FALSE])
      lp_clin[te] <- as.numeric(Xc[te, , drop = FALSE] %*% stats::coef(cfit))
    }
    Cmat[r, "signature"] <- harrell_c(tt, ee, lp_sig)
    Cmat[r, "clinical"]  <- harrell_c(tt, ee, lp_clin)
    Cmat[r, "combined"]  <- harrell_c(tt, ee, lp_comb)
  }

  q <- function(v) stats::quantile(v, c(.025, .975), na.rm = TRUE)
  cindex <- data.frame(
    model = colnames(Cmat),
    cindex = colMeans(Cmat, na.rm = TRUE),
    lower = apply(Cmat, 2L, function(v) q(v)[1L]),
    upper = apply(Cmat, 2L, function(v) q(v)[2L]),
    row.names = NULL
  )
  d <- Cmat[, "combined"] - Cmat[, "clinical"]
  delta <- list(
    mean = mean(d, na.rm = TRUE), lower = unname(q(d)[1L]), upper = unname(q(d)[2L]),
    win_rate = mean(d > 0, na.rm = TRUE)
  )

  structure(list(
    cindex = cindex, delta = delta, per_repeat = Cmat,
    meta = list(n = n, n_folds = n_folds, n_repeats = n_repeats,
      n_predictors = length(predictors), n_clinical = length(clinical), alpha = alpha)
  ), class = "cv_survival_validation")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' @export
print.cv_survival_validation <- function(x, ...) {
  m <- x$meta
  cat(sprintf(
    "Repeated CV survival validation | %d x %d-fold | n=%d | %d genomic vs %d clinical | alpha=%g\n",
    m$n_repeats, m$n_folds, m$n, m$n_predictors, m$n_clinical, m$alpha
  ))
  ci <- x$cindex
  ci[, c("cindex", "lower", "upper")] <- round(ci[, c("cindex", "lower", "upper")], 3)
  cat("\nOut-of-fold C-index (mean [2.5%, 97.5%] across repeats):\n")
  print(ci, row.names = FALSE)
  cat(sprintf(
    "\nAdding genomics to clinical: dC = %+.3f [%+.3f, %+.3f], wins %d%% of repeats\n",
    x$delta$mean, x$delta$lower, x$delta$upper, round(100 * x$delta$win_rate)
  ))
  if (x$delta$lower <= 0) {
    cat("=> interval includes 0: no reliable improvement over clinical alone.\n")
  } else {
    cat("=> interval above 0: genomics adds prognostic value beyond clinical.\n")
  }
  invisible(x)
}
