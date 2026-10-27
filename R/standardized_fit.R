#' @title Fitting quadratic, negative exponential and Jenss-Bayley models on the standardized scale
#' @description Models with these functional forms are badly scaled on raw outcome units: the
#' information matrix of a random-ratio model can have a condition number above 1e8, so the optimizer
#' stops early and the numerical Hessian, and with it the standard errors, depend on rounding. The
#' fitting functions therefore fit these models to standardized outcomes and convert the solution back:
#' \itemize{
#'   \item each outcome process is centred at the mean and divided by the standard deviation of its
#'   first recorded occasion; a time-varying covariate (TVC) is divided by its first-occasion standard
#'   deviation without centring; time is divided by the mean of its last recorded occasion without
#'   centring, so that it runs from 0 to about 1; time-invariant covariates are left as they are, also
#'   when a covariate is one of the time, outcome or TVC columns (for example baseline age as a
#'   growth or class-membership covariate): such a covariate is read from an unscaled copy;
#'   \item the model is built from the same starting values as on the raw scale, mapped exactly to the
#'   standardized scale, and optimized with the same optimizer and retry settings;
#'   \item every free parameter is mapped back with its exact unit factor (and the intercept means
#'   with the centre), the Hessian, its inverse (the parameter covariance matrix), the standard errors
#'   and the gradient with the matching factors, and the -2 log-likelihood is evaluated on the raw data.
#' }
#' The mapping is a reparameterization of the same model, so the maximum-likelihood solution and the
#' likelihood are the same; only the numerical conditioning of the optimization changes. Derived
#' quantities (OpenMx algebras and estimate tables) are computed from the back-transformed model.
#' Set \code{options(nlpsem.standardize = FALSE)} to fit on the raw scale.
#' @keywords internal
#' @noRd
NULL

.std_forms <- c("quadratic", "QUAD", "negative exponential", "EXP", "Jenss-Bayley", "JB")

## Outcome and time exponents of each growth factor (index 0, 1, 2, or g for the ratio).
.std_factor_units <- function(curveFun){
  cf <- c("quadratic" = "QUAD", "QUAD" = "QUAD", "negative exponential" = "EXP", "EXP" = "EXP",
          "Jenss-Bayley" = "JB", "JB" = "JB")[[as.character(curveFun)[1]]]
  list(QUAD = list("0" = c(1, 0), "1" = c(1, -1), "2" = c(1, -2)),
       EXP = list("0" = c(1, 0), "1" = c(1, 0), "g" = c(0, -1)),
       JB = list("0" = c(1, 0), "1" = c(1, -1), "2" = c(1, 0), "g" = c(0, -1)))[[cf]]
}

## Which manifest variables are standardized for a fitting call, with centres and scales; NULL when the
## call is outside the scope (other functional forms, mediation models) or the option is off. Outcomes are
## centred and scaled at their first recorded occasion, a TVC is scaled without centring, and time is
## divided by the mean of its last recorded occasion without centring. A covariate column (growth or
## class-membership TIC) that is also one of these columns gets an unscaled copy under a new name.
.std_plan <- function(dat, sub_Model, curveFun, y_var, records, t_var, TVC = NULL, decompose = NULL,
                      covariates = NULL){
  if (!isTRUE(getOption("nlpsem.standardize", TRUE))) return(NULL)
  if (!isTRUE(as.character(curveFun)[1] %in% .std_forms)) return(NULL)
  if (!isTRUE(sub_Model %in% c("LGCM", "LCSM", "TVC", "MGM"))) return(NULL)
  processes <- list()
  add <- function(key, var, recs, center){
    columns <- paste0(var, recs)
    if (!all(columns %in% names(dat))) return(FALSE)
    first <- suppressWarnings(as.numeric(dat[[columns[1]]]))
    scale <- stats::sd(first, na.rm = TRUE)
    centre <- if (center) mean(first, na.rm = TRUE) else 0
    if (!is.finite(scale) || scale <= 0 || !is.finite(centre)) return(FALSE)
    processes[[key]] <<- list(variable = var, columns = columns, center = centre, scale = scale, centered = center)
    TRUE
  }
  recs <- if (is.list(records)) records else list(records)
  tv <- rep_len(as.character(t_var), length(recs))
  time_columns <- unique(unlist(lapply(seq_along(recs), function(k) paste0(tv[k], recs[[k]]))))
  last <- unique(vapply(seq_along(recs), function(k) paste0(tv[k], utils::tail(recs[[k]], 1)), character(1)))
  if (!all(c(time_columns, last) %in% names(dat))) return(NULL)
  time_scale <- max(vapply(last, function(col) mean(suppressWarnings(as.numeric(dat[[col]])), na.rm = TRUE), numeric(1)))
  if (!is.finite(time_scale) || time_scale <= 0) return(NULL)
  ok <- if (sub_Model == "MGM"){
    all(vapply(seq_along(y_var), function(k) add(y_var[k], y_var[k], records[[k]], TRUE), logical(1)))
  } else {
    add("Y", y_var, records, TRUE) && (sub_Model != "TVC" || add("X", TVC, records, FALSE))
  }
  if (!ok) return(NULL)
  scaled <- c(unlist(lapply(processes, `[[`, "columns"), use.names = FALSE), time_columns)
  shared <- intersect(unique(as.character(covariates)), scaled)
  copies <- make.unique(c(names(dat), paste0("nlpsemTIC_", shared)), sep = "_")
  copies <- stats::setNames(copies[length(dat) + seq_along(shared)], shared)
  list(processes = processes, sub_Model = sub_Model, curveFun = as.character(curveFun)[1], copies = copies,
       decompose = if (is.null(decompose)) 0 else decompose,
       time = list(variables = unique(tv), columns = time_columns, scale = time_scale),
       pairs = if (sub_Model == "MGM" && length(y_var) > 1){
         idx <- utils::combn(seq_along(y_var), 2)
         stats::setNames(lapply(seq_len(ncol(idx)), function(i) y_var[idx[, i]]),
                         apply(idx, 2, function(i) paste0(y_var[i], collapse = "")))
       } else list())
}

## The standardized copy of the data.
.std_data <- function(dat, plan){
  for (col in names(plan$copies)) dat[[plan$copies[[col]]]] <- dat[[col]]
  for (p in plan$processes){
    for (col in p$columns) dat[[col]] <- (dat[[col]] - p$center) / p$scale
  }
  for (col in plan$time$columns) dat[[col]] <- dat[[col]] / plan$time$scale
  dat
}

## Covariate column names as the standardized model reads them (see copies in .std_plan).
.std_cov <- function(plan, columns){
  if (is.null(columns) || !length(plan$copies)) return(columns)
  hit <- columns %in% names(plan$copies)
  columns[hit] <- plan$copies[columns[hit]]
  columns
}

## Exact unit factor (raw = factor * standardized + shift) of every free-parameter label. When a label
## is not recognised, factor is NULL and the labels are listed in unknown, so the call is fitted on the raw scale.
.std_param_map <- function(labels, plan, dimensionless = character()){
  s <- vapply(plan$processes, `[[`, numeric(1), "scale")
  ctr <- vapply(plan$processes, function(p) if (isTRUE(p$centered)) p$center else 0, numeric(1))
  b <- plan$time$scale
  units <- .std_factor_units(plan$curveFun)
  # factor of growth factor `index` of process `key`
  gf <- function(index, key){ u <- units[[index]]; if (is.null(u)) NA_real_ else s[[key]]^u[1] * b^u[2] }
  single <- names(plan$processes)
  prefixes <- c(single, names(plan$pairs))
  prefixes <- prefixes[order(-nchar(prefixes))]
  sx <- if ("X" %in% single) s[["X"]] else NA_real_
  factor <- shift <- stats::setNames(rep(NA_real_, length(labels)), labels)
  for (lab in labels){
    bl <- sub("^c[0-9]+", "", lab)
    f <- NA_real_; sh <- 0
    hit <- prefixes[startsWith(bl, paste0(prefixes, "_"))][1]
    if (lab %in% dimensionless){
      f <- 1
    } else if (identical(hit, "X") && plan$sub_Model == "TVC"){
      rest <- substring(bl, 3L)
      # the TVC process is a nonparametric latent change score model: intercept, slope, relative rates
      xu <- list("0" = c(1, 0), "1" = c(1, -1))
      xf <- function(index){ u <- xu[[index]]; if (is.null(u)) NA_real_ else sx^u[1] * b^u[2] }
      if (grepl("^mueta[01]$", rest)) f <- xf(substr(rest, 6, 6))
      else if (grepl("^psi[01][01]$", rest)) f <- xf(substr(rest, 4, 4)) * xf(substr(rest, 5, 5))
      else if (rest == "residuals") f <- sx^2
      else if (grepl("^rel_rate[0-9]+$", rest)) f <- 1
    } else if (!is.na(hit)){
      rest <- substring(bl, nchar(hit) + 2L)
      if (hit %in% single){
        if (grepl("^mueta[0-9]$", rest)){
          f <- gf(substr(rest, 6, 6), hit); if (rest == "mueta0") sh <- ctr[[hit]]
        } else if (grepl("^psi[0-9g][0-9g]$", rest)){
          f <- gf(substr(rest, 4, 4), hit) * gf(substr(rest, 5, 5), hit)
        } else if (rest == "residuals"){
          f <- s[[hit]]^2
        } else if (rest %in% c("slp_ratio", "acc_ratio")){
          f <- 1 / b
        }
      } else {
        pr <- plan$pairs[[hit]]
        if (grepl("^psi[0-9g][0-9g]$", rest)){
          f <- gf(substr(rest, 4, 4), pr[1]) * gf(substr(rest, 5, 5), pr[2])
        } else if (rest == "RES"){
          f <- s[[pr[1]]] * s[[pr[2]]]
        }
      }
    } else if ("Y" %in% single){
      if (grepl("^beta[0-9g][0-9]+$", bl) || grepl("^beta[0-9g]TIC[0-9]+$", bl)){
        f <- gf(substr(bl, 5, 5), "Y")
      } else if (grepl("^beta[0-9g]TVC$", bl)){
        f <- gf(substr(bl, 5, 5), "Y") / sx
      } else if (bl == "kappa"){
        # effect of the TVC (decompose 0), of its interval-specific slopes (1, rates per unit of time),
        # or of its interval-specific or cumulative changes (2, 3)
        f <- s[["Y"]] / sx * (if (isTRUE(plan$decompose == 1)) b else 1)
      } else if (bl == "Cov_XYres"){
        f <- s[["Y"]] * sx
      } else if (grepl("^covBL[0-9]+$", bl) || grepl("^TVC_m[0-9]+$", bl)){
        f <- sx
      } else if (grepl("^TVC_v[0-9]+$", bl)){
        f <- sx^2
      } else if (grepl("^(mux|phi)[0-9]+$", bl)){
        f <- 1
      }
    }
    factor[[lab]] <- f; shift[[lab]] <- sh
  }
  unknown <- labels[!is.finite(factor)]
  if (length(unknown)) return(list(factor = NULL, shift = NULL, unknown = unknown))
  list(factor = factor, shift = shift, unknown = character())
}

## Free labels of the container of a multiple-group or mixture model (class weights and class-membership
## logit coefficients); all growth-model parameters live in the class submodels. These are dimensionless.
.std_container_labels <- function(model){
  if (!length(model@submodels)) return(character())
  unique(unlist(lapply(model@matrices, function(x){ l <- x@labels[x@free]; l[!is.na(l)] }), use.names = FALSE))
}

## Fit model_mx (built on the raw data) with run(). For calls in the standardized scope (plan not NULL),
## build(plan) returns the same model built on .std_data(dat, plan); it is started from the raw starting
## values mapped exactly to the standardized scale, fitted with the same run(), and converted back.
.fit_model <- function(model_mx, run, plan = NULL, build = NULL){
  if (is.null(plan) || is.null(build) || !methods::is(model_mx, "MxModel")){
    return(list(model = run(model_mx), standardization = NULL))
  }
  start <- OpenMx::omxGetParameters(model_mx)
  if (!length(start) || is.null(names(start)) || anyNA(names(start))){
    return(list(model = run(model_mx), standardization = NULL))
  }
  map <- .std_param_map(names(start), plan, .std_container_labels(model_mx))
  if (is.null(map$factor)){
    return(list(model = run(model_mx),
                standardization = list(applied = FALSE, reason = paste("unrecognised parameter labels:",
                                                                       paste(map$unknown, collapse = ", ")))))
  }
  start_std <- (start - map$shift[names(start)]) / map$factor[names(start)]
  model_std <- OpenMx::omxSetParameters(build(plan), labels = names(start), values = start_std)
  fit_std <- run(model_std)
  model <- .std_back_transform(model_mx, fit_std, map, start, start_std)
  processes <- do.call(rbind, lapply(names(plan$processes), function(k){
    p <- plan$processes[[k]]
    data.frame(process = k, variable = p$variable, first_column = p$columns[1], center = p$center,
               scale = p$scale, centered = p$centered, stringsAsFactors = FALSE)
  }))
  H <- fit_std@output$hessian
  list(model = model,
       standardization = list(applied = TRUE, processes = processes,
                              time = data.frame(variable = plan$time$variables, scale = plan$time$scale,
                                                stringsAsFactors = FALSE),
                              minus2LL_standardized = if (length(fit_std@output)) fit_std@output$Minus2LogLikelihood else NA_real_,
                              condition_number_standardized = .std_condition(H),
                              condition_number_raw_scale = .std_condition(model@output$hessian)))
}

.std_condition <- function(H){
  if (is.null(H) || !all(is.finite(H))) return(NA_real_)
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0) return(NA_real_)
  max(ev) / min(ev)
}

## The raw-scale model at the back-transformed estimates. Its output carries the standardized fit's
## status, the exact transforms of its Hessians, parameter covariance matrix, standard errors and
## gradient, and the -2 log-likelihood evaluated on the raw data.
.std_back_transform <- function(model_mx, fit_std, map, start, start_std){
  est_std <- OpenMx::omxGetParameters(fit_std)
  lab <- names(est_std); f <- map$factor[lab]
  # start + f * (estimate - start) is the exact inverse map and returns the starting values unchanged
  # when the optimizer did not move.
  est <- start[lab] + f * (est_std - start_std[lab])
  model <- OpenMx::omxSetParameters(model_mx, labels = lab, values = est)
  if (!length(fit_std@output)) return(model)
  model <- mxRun(model, useOptimizer = FALSE, silent = TRUE, suppressWarnings = TRUE)
  if (!length(model@output)) return(model)
  so <- fit_std@output; out <- model@output
  scale_vec <- function(x){
    if (is.null(x)) return(NULL)
    nm <- if (is.matrix(x)) rownames(x) else names(x)
    if (is.null(nm) || !all(nm %in% lab)) return(NULL)
    x * f[nm]
  }
  scale_mat <- function(x, power){
    if (is.null(x) || is.null(rownames(x)) || !all(rownames(x) %in% lab) || !identical(rownames(x), colnames(x))) return(NULL)
    x * outer(f[rownames(x)], f[colnames(x)])^power
  }
  keep <- setdiff(names(so), c("estimate", "standardErrors", "hessian", "calculatedHessian", "vcov", "gradient",
                               "Minus2LogLikelihood", "minimum", "fit", "conditionNumber", "matrices", "algebras",
                               "expectations", "data"))
  out[keep] <- so[keep]
  out$estimate <- est[names(so$estimate)]
  out["standardErrors"] <- list(scale_vec(so$standardErrors))
  out["gradient"] <- list(if (is.null(so$gradient)) NULL else so$gradient / f[names(so$gradient)])
  out["hessian"] <- list(scale_mat(so$hessian, -1))
  out["calculatedHessian"] <- list(scale_mat(so$calculatedHessian, -1))
  out["vcov"] <- list(scale_mat(so$vcov, 1))
  out["conditionNumber"] <- list(if (is.null(so$conditionNumber)) NULL else .std_condition(out$hessian))
  model@output <- out
  model
}
