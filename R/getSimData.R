#' @title Simulate Longitudinal Data from the Models of the Package
#'
#' @description This function generates wide-format longitudinal data in which each person's repeated measures
#' follow a latent growth curve, optionally with time-invariant covariates (TICs), a time-varying covariate (TVC),
#' several outcomes, a mediation structure, and several latent classes. The same data can be fitted with a
#' latent growth curve model or a latent change score model (\code{getLGCM()} or \code{getLCSM()}, or
#' \code{y_model} in \code{getTVCmodel()}, \code{getMGM()}, \code{getMGroup()} and \code{getMIX()}), and
#' mediation data with \code{getMediation()}. The true parameter values are attached to the data under the names
#' that these fitting functions use in their \code{Estimates}, which are the same for both models, so that a fit
#' can be compared with the values that generated the data.
#'
#' @param n An integer: the number of persons. Default is \code{500}.
#' @param sub_Model A string: \code{"LGCM"} (one outcome), \code{"TVC"} (an outcome with a time-varying
#' covariate), \code{"MGM"} (several outcomes) or \code{"MED"} (longitudinal mediation). Default is \code{"LGCM"}.
#' @param curveFun A string specifying the functional form of the growth curves: \code{"linear"} (or
#' \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or \code{"EXP"}),
#' \code{"Jenss-Bayley"} (or \code{"JB"}), \code{"bilinear spline"} (or \code{"BLS"}) or \code{"nonparametric"}
#' (or \code{"NonP"}); for \code{sub_Model = "MED"}, \code{"linear"} or \code{"bilinear spline"}. A latent growth
#' curve model fits all but the nonparametric form, and a latent change score model fits the quadratic, negative
#' exponential, Jenss-Bayley and nonparametric forms.
#' @param intrinsic A logical flag: whether the growth ratio (negative exponential and Jenss-Bayley) or the knot
#' (bilinear spline) varies between persons. Not used for \code{sub_Model = "MED"}, whose knots are the same for
#' everyone. Default is \code{FALSE}.
#' @param nT An integer: the number of measurement occasions. Default is \code{10}.
#' @param times A numeric vector of the \code{nT} planned occasions. Default is \code{0, 1, ..., nT - 1}.
#' @param jitter A non-negative number: each person's interior occasions are drawn uniformly within
#' \code{jitter} of the planned ones (individually varying measurement occasions); the first and the last
#' occasions are the planned ones. Default is \code{0.25}.
#' @param decompose For \code{sub_Model = "TVC"}: \code{0} (the TVC enters directly), \code{1} (through its
#' interval-specific slopes), \code{2} (its interval-specific changes) or \code{3} (its change from baseline).
#' @param x_type For \code{sub_Model = "MED"}: \code{"baseline"} or \code{"longitudinal"} predictor.
#' @param params A list of parameters for one class, or a list of such lists, one per latent class. See
#' Details. Default is \code{NULL}: one class with the settings of the package's simulation studies.
#' @param prop A numeric vector of class proportions for several classes without cluster TICs. Default is
#' \code{NULL} (equal proportions).
#' @param cluster A list for several classes whose membership depends on cluster TICs: \code{mean} and
#' \code{cov} of the cluster TICs, and \code{beta}, a matrix of multinomial-logit coefficients with one row per
#' class (the first row zeros, for the reference class) and one column for the intercept and each TIC.
#' Default is \code{NULL}.
#' @param growth_TIC A character vector of column names for the growth TICs (\code{"LGCM"} and \code{"TVC"}). Default is \code{"ex1"}, \code{"ex2"}, ... when \code{params} has growth-TIC effects.
#' @param cluster_TIC A character vector of column names for the cluster TICs. Default is \code{"gx1"},
#' \code{"gx2"}, ... when \code{cluster} is given.
#' @param y_var The prefix of the outcome columns; for \code{sub_Model = "MGM"}, one prefix per process, with the
#' restrictions of \code{\link{getMGM}}. Defaults are \code{"Y"} and, for \code{"MGM"}, \code{c("Y", "Z")}.
#' @param t_var The prefix of the occasion columns, shared by all processes. Default is \code{"T"}.
#' @param TVC The prefix of the TVC columns. Default is \code{"TVC"}.
#' @param m_var,x_var The prefixes of the mediator and predictor columns for \code{sub_Model = "MED"} (a baseline
#' predictor is the single column \code{x_var}). Defaults are \code{"M"} and \code{"X"}.
#' @param seed An integer seed for reproducible data, or \code{NULL}. The random number generator state of the
#' session is restored afterwards. Default is \code{NULL}.
#'
#' @return A data frame with columns \code{id}, the repeated measures, the occasions \code{t_var1}, ..., the
#' TICs and, for several classes, \code{class}. Its attribute \code{"truth"} holds the true parameter values,
#' named as in the \code{Estimates} of the matching fitting function, or of \code{getMIX()} for several classes
#' (with the class weights as log odds \code{p2}, ... against the first class, or as cluster-TIC coefficients
#' \code{beta20}, ...); derived quantities (rates and changes of a latent change score model, mediation and
#' total effects) are not included. Attribute \code{"growth_factors"} holds each person's growth factors, and
#' attribute \code{"settings"} the arguments.
#'
#' @details Growth factors, TICs and the baseline of a TVC or of a mediation predictor are drawn from
#' multivariate normal distributions, each person's repeated measures follow the functional form at their own
#' occasions, and normal residuals are added (correlated between processes at the same occasion). Each
#' trajectory is the exact growth curve at the person's own growth factors; a latent change score model
#' describes the same trajectory by its rates of change. The nonparametric form is a growth curve whose slope in
#' each interval is \code{eta1} times the interval's relative rate.
#'
#' Two approximations in the fitting functions can affect recovery of the generating values. A model with a
#' random growth ratio or knot uses a first-order Taylor expansion around its mean. A latent change score model
#' approximates each interval's change by its midpoint rate times its length. The latter is exact for quadratic
#' and nonparametric curves, but not for negative exponential and Jenss-Bayley curves; its error shrinks with
#' shorter intervals or less curvature.
#'
#' In one generated data set of 5,000 persons per configuration, using the default parameters and ten yearly
#' occasions, the observed errors in growth-factor means were about 1\% for models using the Taylor expansion
#' alone. For a negative exponential TVC model with \code{decompose = 1} or \code{3}, the TVC effect was
#' 5\% to 10\% above its generating value and the mean of \code{eta1} or of the ratio up to 3\% below.
#' The latent change score fits had errors up to about 2.5\% in the magnitude of the mean of \code{eta1}
#' (negative exponential) or \code{eta2} (Jenss-Bayley); the random-ratio Jenss-Bayley fit's covariance of
#' \code{eta2} and the ratio was about 40\% above its generating value. Matched-data comparisons using the
#' fitting approximations reduced these discrepancies. These percentages describe those particular samples,
#' not a replicated estimate of bias, a bound on error, or an assessment of interval coverage. The reproducible
#' examples are in \file{data-raw/approximation-bias} in the source repository.
#'
#' The elements of a class list of \code{params} are:
#' \describe{
#'   \item{\code{"LGCM"}}{\code{mean}: a named numeric vector of the means of the growth factors:
#'   \code{eta0} (intercept), \code{eta1} (slope; for the negative exponential, the change to the asymptote; for
#'   the bilinear spline, the slope before the knot; for the nonparametric form, the slope of the first
#'   interval), \code{eta2} (quadratic coefficient; for the Jenss-Bayley, the coefficient of the exponential
#'   term; for the bilinear spline, the slope after the knot) and \code{gamma} (growth ratio or knot), as the
#'   form requires; \code{cov}: the covariance matrix of the growth factors that vary between persons (in the
#'   order above), not explained by the growth TICs; \code{residual}: the residual variance; \code{rel_rate}:
#'   for the nonparametric form, the relative rates of the \code{nT - 1} intervals, the first \code{1};
#'   \code{tic_mean}, \code{tic_cov}, \code{beta}: for growth TICs, their means and covariance matrix, and their
#'   effects on the growth factors that vary between persons (one row per growth factor, one column per TIC).}
#'   \item{\code{"MGM"}}{\code{mean}: a list of mean vectors as above, one per process; \code{cov}: the
#'   covariance matrix of the growth factors that vary between persons, process by process; \code{residual}: the
#'   covariance matrix of the residuals of the processes at the same occasion; \code{rel_rate}: for the
#'   nonparametric form, a list of relative rates, one per process.}
#'   \item{\code{"TVC"}}{For the outcome, the elements of a latent growth curve or change score model as above,
#'   and \code{kappa}, the effect of the TVC. With \code{decompose = 0}, \code{tvc}: a list with the \code{mean}
#'   and \code{var} of the TVC at each occasion. With \code{decompose = 1}, \code{2} or \code{3}, the TVC is a
#'   nonparametric latent change score process: \code{tvc} is a list with its \code{mean} (\code{eta0},
#'   \code{eta1}), \code{cov}, \code{rel_rate} and \code{residual}; \code{beta_tvc}: the effects of the TVC's
#'   latent baseline on the outcome's growth factors that vary between persons; \code{res_cov}: the covariance of
#'   the TVC and outcome residuals at the same occasion; with growth TICs, \code{cov_bl}: their covariances
#'   with the TVC's latent baseline. The outcome's \code{mean} holds marginal means.}
#'   \item{\code{"MED"}}{\code{x}: the predictor, \code{list(mean, var)} at baseline or a growth process
#'   \code{list(mean, cov)}; \code{m}, \code{y}: the mediator and outcome processes, \code{list(mean, cov)}
#'   (marginal means, and covariances not explained by the predictor and mediator); \code{beta_xm},
#'   \code{beta_xy}: the effects of the predictor (or of its growth factors) on the mediator's and outcome's
#'   growth factors, and \code{beta_my}: of the mediator's growth factors on the outcome's (one row per
#'   receiving growth factor; a growth factor affects only itself and later ones); \code{residual}: the
#'   covariance matrix of the residuals of the (predictor,) mediator and outcome at the same occasion. The growth
#'   factors are \code{eta0} and \code{eta1} for the linear form, and \code{eta1} (slope before the knot),
#'   \code{etar} (status at the knot), \code{eta2} (slope after the knot) and the knot \code{gamma} for the
#'   bilinear spline.}
#' }
#'
#' @examples
#' dat <- getSimData(n = 200, curveFun = "QUAD", seed = 20261027)
#' head(dat[, 1:6])
#' attr(dat, "truth")
#' \donttest{
#' fit <- getLGCM(dat = dat, t_var = "T", y_var = "Y", curveFun = "QUAD", intrinsic = FALSE,
#'                records = 1:10, paramOut = TRUE)
#' truth <- attr(dat, "truth")
#' cbind(estimate = fit@Estimates$Estimate[match(names(truth), fit@Estimates$Name)], truth)
#' }
#' bivariate <- getSimData(n = 200, sub_Model = "MGM", curveFun = "LIN", seed = 1)
#' names(attr(bivariate, "truth"))
#'
#' @importFrom stats rnorm runif
#' @export
getSimData <- function(n = 500, sub_Model = "LGCM", curveFun, intrinsic = FALSE, nT = 10, times = NULL,
                       jitter = 0.25, decompose = NULL, x_type = NULL, params = NULL, prop = NULL,
                       cluster = NULL, growth_TIC = NULL, cluster_TIC = NULL, y_var = NULL, t_var = "T",
                       TVC = "TVC", m_var = "M", x_var = "X", seed = NULL){
  validate_curveFun(curveFun)
  if (!isTRUE(sub_Model %in% c("LGCM", "TVC", "MGM", "MED"))){
    stop("sub_Model must be 'LGCM', 'TVC', 'MGM' or 'MED'; the data follow growth curves, which a latent change ",
         "score model can also fit.")
  }
  form <- .sim_form(curveFun)
  # Growth curves throughout; only the nonparametric form takes the names of a latent change score model.
  y_model <- if (form == "NonP") "LCSM" else "LGCM"
  if (sub_Model == "MED"){
    if (!form %in% c("LIN", "BLS")) stop("A mediation model has a linear or bilinear spline functional form.")
    if (!isTRUE(x_type %in% c("baseline", "longitudinal"))) stop("x_type must be 'baseline' or 'longitudinal'.")
    intrinsic <- FALSE
  } else{
    validate_curveFun_ymodel(curveFun, y_model)
    validate_intrinsic(intrinsic, curveFun, model_type = y_model)
  }
  if (sub_Model == "TVC" && !isTRUE(decompose %in% 0:3)) stop("decompose must be 0, 1, 2 or 3.")
  if (is.null(y_var)) y_var <- if (sub_Model == "MGM") c("Y", "Z") else "Y"
  if (sub_Model == "MGM" && (length(y_var) < 2L || anyDuplicated(y_var))) stop("'y_var' must name two or more processes.")
  if (sub_Model == "MGM") validate_mgm_names(y_var, y_model)
  if (sub_Model != "MGM" && length(y_var) != 1L) stop("'y_var' must be one prefix.")
  for (nm in c("t_var", "TVC", "m_var", "x_var")){
    value <- get(nm)
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)){
      stop("'", nm, "' must be one nonempty character prefix.")
    }
    assign(nm, unname(value))
  }
  if (!is.character(y_var) || anyNA(y_var) || any(!nzchar(y_var))) stop("'y_var' must be nonempty character prefixes.")
  y_var <- unname(y_var)
  check_tic_names <- function(x, what){
    if (!is.null(x) && (!is.character(x) || anyNA(x) || any(!nzchar(x)))){
      stop("'", what, "' must contain nonempty character column names.")
    }
    unname(x)
  }
  growth_TIC <- check_tic_names(growth_TIC, "growth_TIC")
  cluster_TIC <- check_tic_names(cluster_TIC, "cluster_TIC")
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 2 || n != round(n)) stop("'n' must be a finite integer of at least 2.")
  if (!is.numeric(nT) || length(nT) != 1L || !is.finite(nT) || nT < 3 || nT != round(nT)) stop("'nT' must be a finite integer of at least 3.")
  if (is.null(times)) times <- seq_len(nT) - 1
  if (!is.numeric(times) || length(times) != nT || any(!is.finite(times)) || any(diff(times) <= 0)){
    stop("'times' must be ", nT, " finite increasing numbers.")
  }
  if (!is.numeric(jitter) || length(jitter) != 1L || !is.finite(jitter) || jitter < 0 || (nT > 2 && jitter >= min(diff(times)) / 2)){
    stop("'jitter' must be non-negative and smaller than half the shortest interval between planned occasions.")
  }
  if (!is.null(seed)){
    if (exists(".Random.seed", envir = globalenv())){
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    } else{
      on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }

  intrinsic <- isTRUE(intrinsic) && form %in% c("EXP", "JB", "BLS")
  ctx <- list(sub_Model = sub_Model, curveFun = curveFun, form = form, intrinsic = intrinsic, y_model = y_model,
              decompose = decompose, x_type = x_type, nT = nT, y_var = y_var)
  ctx$factors <- if (sub_Model == "MED") .sim_med_factors(form) else .sim_factors(form)
  ctx$random <- if (sub_Model == "MED") setdiff(ctx$factors, "gamma") else
    if (intrinsic || !"gamma" %in% ctx$factors) ctx$factors else setdiff(ctx$factors, "gamma")
  if (is.null(params)) params <- list(.sim_default_params(ctx))
  if (!.sim_is_class_list(params, sub_Model)) params <- list(params)
  nClass <- length(params)
  params <- lapply(seq_len(nClass), function(k) .sim_check_params(params[[k]], k, ctx))
  nTIC <- unique(vapply(params, function(p) if (is.null(p[["beta"]])) 0L else ncol(p[["beta"]]), integer(1)))
  if (length(nTIC) != 1L) stop("All classes need the same number of growth TICs.")
  if (nTIC > 0 && is.null(growth_TIC)) growth_TIC <- paste0("ex", seq_len(nTIC))
  if (length(growth_TIC) != nTIC) stop("'growth_TIC' must name the ", nTIC, " growth TIC(s) in 'params'.")
  ctx$growth_TIC <- growth_TIC

  # Class membership
  X <- NULL
  if (nClass == 1L){
    if (!is.null(prop) || !is.null(cluster)) stop("'prop' and 'cluster' need parameters for several classes.")
    class <- rep(1L, n)
  } else if (!is.null(cluster)){
    if (!is.null(prop)) stop("Give either 'prop' or 'cluster', not both.")
    nX <- length(cluster$mean)
    if (!is.numeric(cluster$mean) || !nX || !all(is.finite(cluster$mean))) stop("'cluster$mean' must be finite numbers.")
    cluster$cov <- as.matrix(cluster$cov)
    if (!identical(dim(cluster$cov), c(nX, nX)) || !is.numeric(cluster$cov) || !all(is.finite(cluster$cov)) ||
        !isSymmetric(unname(cluster$cov)) || min(eigen(cluster$cov, symmetric = TRUE, only.values = TRUE)$values) <= 0){
      stop("'cluster$cov' must be a symmetric positive-definite ", nX, " x ", nX, " matrix.")
    }
    if (!is.matrix(cluster$beta) || !identical(dim(cluster$beta), c(nClass, nX + 1L)) || !all(is.finite(cluster$beta)) ||
        any(cluster$beta[1, ] != 0)){
      stop("'cluster$beta' must be a ", nClass, " x ", nX + 1L, " matrix of finite numbers with a first row of zeros.")
    }
    if (is.null(cluster_TIC)) cluster_TIC <- paste0("gx", seq_len(nX))
    if (length(cluster_TIC) != nX) stop("'cluster_TIC' must name the ", nX, " cluster TIC(s).")
    X <- .sim_mvn(n, cluster$mean, cluster$cov, "cluster$cov")
    eta <- cbind(1, X) %*% t(cluster$beta)
    prob <- exp(eta - apply(eta, 1, max)); prob <- prob / rowSums(prob)
    class <- 1L + rowSums(runif(n) > t(apply(prob, 1, cumsum))[, -nClass, drop = FALSE])
  } else{
    if (is.null(prop)) prop <- rep(1 / nClass, nClass)
    if (!is.numeric(prop) || length(prop) != nClass || !all(is.finite(prop)) || any(prop <= 0)) stop("'prop' must be ", nClass, " positive numbers.")
    prop <- prop / sum(prop)
    class <- sample.int(nClass, n, replace = TRUE, prob = prop)
  }

  # Column names of the data: the processes are stored under their roles and named only here
  roles <- if (sub_Model == "MED") c(if (x_type == "longitudinal") c(.X = x_var), .M = m_var, .y1 = y_var) else
    c(stats::setNames(y_var, paste0(".y", seq_along(y_var))), if (sub_Model == "TVC") c(.TVC = TVC))
  out_names <- c("id", unlist(lapply(roles, paste0, seq_len(nT)), use.names = FALSE), paste0(t_var, seq_len(nT)),
                 if (sub_Model == "MED" && x_type == "baseline") x_var, growth_TIC, if (!is.null(X)) cluster_TIC,
                 if (nClass > 1L) "class")
  if (anyDuplicated(out_names)){
    stop("The prefixes and covariate names give duplicated columns: ",
         paste(unique(out_names[duplicated(out_names)]), collapse = ", "), ".")
  }

  # Individual occasions
  occasions <- matrix(times, n, nT, byrow = TRUE)
  if (jitter > 0) occasions[, 2:(nT - 1)] <- occasions[, 2:(nT - 1)] + runif(n * (nT - 2), -jitter, jitter)

  # Repeated measures, TICs and growth factors, class by class
  generate <- switch(EXPR = sub_Model, LGCM = .sim_class_uni, TVC = .sim_class_tvc, MGM = .sim_class_mgm,
                     MED = .sim_class_med)
  columns <- list(); GF <- NULL
  for (k in seq_len(nClass)){
    rows <- which(class == k)
    if (!length(rows)) next
    out <- generate(params[[k]], occasions[rows, , drop = FALSE], ctx, k)
    for (nm in names(out$columns)){
      if (is.null(columns[[nm]])) columns[[nm]] <- matrix(NA_real_, n, ncol(out$columns[[nm]]))
      columns[[nm]][rows, ] <- out$columns[[nm]]
    }
    if (is.null(GF)) GF <- matrix(NA_real_, n, ncol(out$gf), dimnames = list(NULL, colnames(out$gf)))
    GF[rows, ] <- out$gf
  }

  dat <- data.frame(id = seq_len(n))
  for (role in names(roles)){
    for (j in seq_len(nT)) dat[[paste0(roles[[role]], j)]] <- columns[[role]][, j]
  }
  for (j in seq_len(nT)) dat[[paste0(t_var, j)]] <- occasions[, j]
  if (!is.null(columns$.XB)) dat[[x_var]] <- columns$.XB[, 1]
  if (nTIC > 0) for (j in seq_len(nTIC)) dat[[growth_TIC[j]]] <- columns$.TIC[, j]
  if (!is.null(X)) for (j in seq_along(cluster_TIC)) dat[[cluster_TIC[j]]] <- X[, j]
  if (nClass > 1L) dat$class <- class

  attr(dat, "truth") <- .sim_truth(ctx, params, prop, cluster)
  attr(dat, "growth_factors") <- GF
  attr(dat, "settings") <- list(n = n, sub_Model = sub_Model, curveFun = curveFun, intrinsic = intrinsic, nT = nT,
                                times = times, jitter = jitter, decompose = decompose,
                                x_type = x_type, params = params, prop = prop, cluster = cluster,
                                growth_TIC = growth_TIC, cluster_TIC = cluster_TIC, y_var = y_var, t_var = t_var,
                                TVC = TVC, m_var = m_var, x_var = x_var, seed = seed)
  dat
}

## Canonical short name of a functional form.
.sim_form <- function(curveFun){
  c(linear = "LIN", LIN = "LIN", quadratic = "QUAD", QUAD = "QUAD", "negative exponential" = "EXP", EXP = "EXP",
    "Jenss-Bayley" = "JB", JB = "JB", "bilinear spline" = "BLS", BLS = "BLS", nonparametric = "NonP",
    NonP = "NonP")[[curveFun]]
}

## Growth factors of a form, in the order of the fitting functions (the ratio or knot last).
.sim_factors <- function(form){
  list(LIN = c("eta0", "eta1"), QUAD = c("eta0", "eta1", "eta2"), EXP = c("eta0", "eta1", "gamma"),
       JB = c("eta0", "eta1", "eta2", "gamma"), BLS = c("eta0", "eta1", "eta2", "gamma"),
       NonP = c("eta0", "eta1"))[[form]]
}

## Growth factors of a mediation model: the bilinear spline is parameterized at the knot.
.sim_med_factors <- function(form){
  if (form == "LIN") c("eta0", "eta1") else c("eta1", "etar", "eta2", "gamma")
}

## Labels of the growth factors in the names of the estimates.
.sim_label <- function(f) c(eta0 = "0", eta1 = "1", eta2 = "2", gamma = "g", etar = "r")[f]

## Whether params is already a list of class lists.
.sim_is_class_list <- function(params, sub_Model){
  key <- if (sub_Model == "MED") "m" else "mean"
  is.list(params) && is.null(params[[key]]) && all(vapply(params, function(p) is.list(p) && !is.null(p[[key]]), logical(1)))
}

## A covariance matrix with the given standard deviations and correlation 0.3.
.sim_cov03 <- function(sd){
  R <- matrix(0.3, length(sd), length(sd)); diag(R) <- 1
  R * outer(sd, sd)
}

## Settings of the package's simulation studies: intercept variance 25, other variances 1 (the growth ratio
## 0.1^2, the knot 0.3^2), correlations 0.3 and residual variance 1.
.sim_default_params <- function(ctx){
  uni <- function(form, intrinsic, nT){
    d <- list(LIN = list(mean = c(eta0 = 50, eta1 = 5), sd = c(5, 1)),
              QUAD = list(mean = c(eta0 = 50, eta1 = 20, eta2 = -1), sd = c(5, 1, 0.3)),
              EXP = list(mean = c(eta0 = 50, eta1 = 30, gamma = 0.6), sd = c(5, 3, 0.1)),
              JB = list(mean = c(eta0 = 50, eta1 = 2.5, eta2 = -30, gamma = -0.7), sd = c(5, 1, 3, 0.1)),
              BLS = list(mean = c(eta0 = 100, eta1 = -5, eta2 = -1.8, gamma = 4.5), sd = c(5, 1, 1, 0.3)),
              NonP = list(mean = c(eta0 = 50, eta1 = 5), sd = c(5, 1)))[[form]]
    sd <- d$sd
    if (!intrinsic && form %in% c("EXP", "JB", "BLS")) sd <- sd[-length(sd)]
    out <- list(mean = d$mean, cov = .sim_cov03(sd), residual = 1, sd = sd)
    if (form == "NonP") out$rel_rate <- .sim_default_rates(nT)
    out
  }
  form <- ctx$form; nT <- ctx$nT
  if (ctx$sub_Model == "LGCM"){
    p <- uni(form, ctx$intrinsic, nT); p$sd <- NULL
    return(p)
  }
  if (ctx$sub_Model == "MGM"){
    base <- uni(form, ctx$intrinsic, nT); nP <- length(ctx$y_var); sd <- rep(base$sd, nP)
    cov <- .sim_cov03(sd)
    nr <- length(base$sd)
    for (q in seq_len(nP)){ i <- (q - 1) * nr + seq_len(nr); cov[i, i] <- base$cov }
    residual <- matrix(0.3, nP, nP); diag(residual) <- 1
    mean <- lapply(seq_len(nP), function(q){ m <- base$mean; m["eta0"] <- m["eta0"] + 4 * (q - 1); m })
    names(mean) <- ctx$y_var
    out <- list(mean = mean, cov = cov, residual = residual)
    if (form == "NonP") out$rel_rate <- stats::setNames(rep(list(base$rel_rate), nP), ctx$y_var)
    return(out)
  }
  if (ctx$sub_Model == "TVC"){
    p <- uni(form, ctx$intrinsic, nT); sd <- p$sd; p$sd <- NULL
    p$kappa <- 0.5
    rates <- .sim_default_rates(nT)
    if (ctx$decompose == 0){
      p$tvc <- list(mean = 5 * c(0, cumsum(rates)), var = rep(4, nT))
    } else{
      p$tvc <- list(mean = c(eta0 = 0, eta1 = 5), cov = .sim_cov03(c(4, 1)), rel_rate = rates, residual = 1)
      p$beta_tvc <- 0.25 * sd / 4
      p$res_cov <- 0.3
    }
    return(p)
  }
  # Mediation
  if (form == "LIN"){
    if (ctx$x_type == "baseline"){
      list(x = list(mean = 0, var = 1),
           m = list(mean = c(eta0 = 50, eta1 = 5), cov = .sim_cov03(c(5, 1))),
           y = list(mean = c(eta0 = 35, eta1 = 5), cov = .sim_cov03(c(5, 1))),
           beta_xm = c(2.55, 0.51), beta_xy = c(2.55, 0.51) / sqrt(2),
           beta_my = matrix(c(3 / 7, 0, 2 / 40, 3 / 7), 2, byrow = TRUE) / sqrt(2),
           residual = matrix(c(1, 0.3, 0.3, 1), 2))
    } else{
      b <- matrix(c(2.55 / 4, 0, 0.102 / 2, 0.51 / 4), 2, byrow = TRUE)
      list(x = list(mean = c(eta0 = 60, eta1 = 5), cov = .sim_cov03(c(5, 1))),
           m = list(mean = c(eta0 = 80, eta1 = 5), cov = .sim_cov03(c(5, 1))),
           y = list(mean = c(eta0 = 100, eta1 = 5), cov = .sim_cov03(c(5, 1))),
           beta_xm = b, beta_xy = b, beta_my = matrix(c(1 / 9, 0, 1 / 18, 1 / 9), 2, byrow = TRUE) * sqrt(2),
           residual = .sim_cov03(c(1, 1, 1)))
    }
  } else{
    my <- matrix(c(3 / 7, 0, 0, 10 / 8, 3 / 7, 0, 2 / 8, 2 / 40, 3 / 7), 3, byrow = TRUE) / sqrt(2)
    if (ctx$x_type == "baseline"){
      list(x = list(mean = 0, var = 1),
           m = list(mean = c(eta1 = 5, etar = 50, eta2 = 2.6, gamma = 3.5), cov = .sim_cov03(c(1, 5, 1))),
           y = list(mean = c(eta1 = 5, etar = 35, eta2 = 2.6, gamma = 5.5), cov = .sim_cov03(c(1, 5, 1))),
           beta_xm = c(0.51, 2.55, 0.51), beta_xy = c(0.51, 2.55, 0.51) / sqrt(2), beta_my = my,
           residual = matrix(c(1, 0.3, 0.3, 1), 2))
    } else{
      b <- matrix(c(0.51 / 4, 0, 0, 2.55 / 4, 0.51 / 4, 0, 0.102 / 2, 0.102 / 2, 0.51 / 4), 3, byrow = TRUE)
      list(x = list(mean = c(eta1 = 5, etar = 60, eta2 = 2.6, gamma = 4.5), cov = .sim_cov03(c(1, 5, 1))),
           m = list(mean = c(eta1 = 5, etar = 80, eta2 = 2.6, gamma = 4.5), cov = .sim_cov03(c(1, 5, 1))),
           y = list(mean = c(eta1 = 5, etar = 100, eta2 = 2.6, gamma = 5.5), cov = .sim_cov03(c(1, 5, 1))),
           beta_xm = b, beta_xy = b,
           beta_my = matrix(c(1 / 9, 0, 0, 1 / 18, 1 / 9, 0, 1 / 18, 1 / 18, 1 / 9), 3, byrow = TRUE) * sqrt(2),
           residual = .sim_cov03(c(1, 1, 1)))
    }
  }
}

## Relative rates of the nT - 1 intervals: from 1 down to 0.2.
.sim_default_rates <- function(nT) seq(1, by = -0.8 / (nT - 2), length.out = nT - 1)

.sim_check_params <- function(p, k, ctx){
  where <- paste0(" (class ", k, ")")
  if (!is.list(p)) stop("Each class needs a list of parameters", where, ".")
  # Values are checked before any coercion: finite numbers of the exact length, symmetric positive-definite
  # covariance matrices, and named means reordered to the order of the growth factors.
  check_num <- function(x, what, len = NULL){
    if (!is.numeric(x) || !length(x) || !all(is.finite(x)) || (!is.null(len) && length(x) != len)){
      stop("'", what, "' must be ", if (is.null(len)) "finite numbers" else paste(len, "finite number(s)"), where, ".")
    }
    x
  }
  check_positive <- function(x, what){
    if (check_num(x, what, 1L) <= 0) stop("'", what, "' must be positive", where, ".")
    x
  }
  check_cov <- function(cov, size, what, labels){
    cov <- as.matrix(cov)
    if (!identical(dim(cov), c(size, size))) stop("'", what, "' must be a ", size, " x ", size, " matrix for ", labels, where, ".")
    check_num(cov, what)
    if (!isSymmetric(unname(cov))) stop("'", what, "' must be symmetric", where, ".")
    if (min(eigen(cov, symmetric = TRUE, only.values = TRUE)$values) <= 0) stop("'", what, "' must be positive definite", where, ".")
    cov
  }
  check_mean <- function(mean, factors, what){
    nm <- names(mean)
    if (!is.numeric(mean) || is.null(nm) || anyDuplicated(nm) || !setequal(nm, factors)){
      stop("'", what, "' must name ", paste(factors, collapse = ", "), " once each", where, ".")
    }
    check_num(mean[factors], what)
  }
  check_rates <- function(rel_rate, what){
    if (is.null(rel_rate)) stop("'", what, "' is missing", where, ".")
    check_num(rel_rate, what, ctx$nT - 1L)
    if (rel_rate[1] != 1) stop("'", what, "' must give the relative rates of the ", ctx$nT - 1L, " intervals, the first 1", where, ".")
    rel_rate
  }
  check_tic <- function(p, nr){
    if (is.null(p[["beta"]])) return(p)
    p[["beta"]] <- as.matrix(check_num(p[["beta"]], "beta"))
    nTIC <- ncol(p[["beta"]])
    if (nrow(p[["beta"]]) != nr || is.null(p$tic_mean) || is.null(p$tic_cov)){
      stop("Growth TICs need 'tic_mean', 'tic_cov' and a 'beta' with one row per growth factor that varies",
           " between persons", where, ".")
    }
    p$tic_mean <- check_num(p$tic_mean, "tic_mean", nTIC)
    p$tic_cov <- check_cov(p$tic_cov, nTIC, "tic_cov", "the growth TICs")
    p
  }
  factors <- ctx$factors; random <- ctx$random; nr <- length(random)
  if (ctx$sub_Model %in% c("LGCM", "TVC")){
    if (is.null(p$mean) || is.null(p$cov) || is.null(p$residual)) stop("Each class needs 'mean', 'cov' and 'residual'", where, ".")
    p$mean <- check_mean(p$mean, factors, "mean")
    p$cov <- check_cov(p$cov, nr, "cov", paste(random, collapse = ", "))
    p$residual <- check_positive(p$residual, "residual")
    if (ctx$form == "NonP") p$rel_rate <- check_rates(p$rel_rate, "rel_rate")
    p <- check_tic(p, nr)
  }
  if (ctx$sub_Model == "TVC"){
    p$kappa <- check_num(p$kappa, "kappa", 1L)
    if (!is.list(p$tvc)) stop("'tvc' must be a list", where, ".")
    if (ctx$decompose == 0){
      p$tvc$mean <- check_num(p$tvc$mean, "tvc$mean", ctx$nT)
      p$tvc$var <- check_num(p$tvc$var, "tvc$var", ctx$nT)
      if (any(p$tvc$var <= 0)) stop("'tvc$var' must be positive at each of the ", ctx$nT, " occasions", where, ".")
    } else{
      p$tvc$mean <- check_mean(p$tvc$mean, c("eta0", "eta1"), "tvc$mean")
      p$tvc$cov <- check_cov(p$tvc$cov, 2L, "tvc$cov", "eta0, eta1")
      p$tvc$rel_rate <- check_rates(p$tvc$rel_rate, "tvc$rel_rate")
      p$tvc$residual <- check_positive(p$tvc$residual, "tvc$residual")
      p$beta_tvc <- check_num(p$beta_tvc, "beta_tvc", nr)
      p$res_cov <- check_num(p$res_cov, "res_cov", 1L)
      nTIC <- if (is.null(p[["beta"]])) 0L else ncol(p[["beta"]])
      if (is.null(p$cov_bl)) p$cov_bl <- rep(0, nTIC)
      if (nTIC > 0) p$cov_bl <- check_num(p$cov_bl, "cov_bl", nTIC) else if (length(p$cov_bl)) stop("'cov_bl' needs growth TICs", where, ".")
    }
  }
  if (ctx$sub_Model == "MGM"){
    nP <- length(ctx$y_var)
    if (!is.list(p$mean) || length(p$mean) != nP) stop("'mean' must be a list with one mean vector per process", where, ".")
    if (is.null(names(p$mean))) names(p$mean) <- ctx$y_var
    if (!setequal(names(p$mean), ctx$y_var) || anyDuplicated(names(p$mean))) stop("'mean' must name the processes ", paste(ctx$y_var, collapse = ", "), where, ".")
    p$mean <- lapply(p$mean[ctx$y_var], check_mean, factors = factors, what = "mean")
    p$cov <- check_cov(p$cov, nP * nr, "cov", paste0(nr, " growth factors of each of ", nP, " processes"))
    if (is.null(dim(p$residual))) p$residual <- diag(check_num(p$residual, "residual", nP), nP)
    p$residual <- check_cov(p$residual, nP, "residual", "the processes")
    if (ctx$form == "NonP"){
      if (!is.list(p$rel_rate) || length(p$rel_rate) != nP) stop("'rel_rate' must be a list with one vector per process", where, ".")
      if (is.null(names(p$rel_rate))) names(p$rel_rate) <- ctx$y_var
      if (!setequal(names(p$rel_rate), ctx$y_var)) stop("'rel_rate' must name the processes ", paste(ctx$y_var, collapse = ", "), where, ".")
      p$rel_rate <- lapply(p$rel_rate[ctx$y_var], check_rates, what = "rel_rate")
    }
  }
  if (ctx$sub_Model == "MED"){
    nf <- length(random)
    long <- ctx$x_type == "longitudinal"
    if (!is.list(p$x) || !is.list(p$m) || !is.list(p$y)) stop("A mediation model needs 'x', 'm' and 'y'", where, ".")
    if (long){
      p$x$mean <- check_mean(p$x$mean, factors, "x$mean")
      p$x$cov <- check_cov(p$x$cov, nf, "x$cov", paste(random, collapse = ", "))
    } else{
      p$x$mean <- check_num(p$x$mean, "x$mean", 1L)
      p$x$var <- check_positive(p$x$var, "x$var")
    }
    for (v in c("m", "y")){
      p[[v]]$mean <- check_mean(p[[v]]$mean, factors, paste0(v, "$mean"))
      p[[v]]$cov <- check_cov(p[[v]]$cov, nf, paste0(v, "$cov"), paste(random, collapse = ", "))
    }
    path <- function(b, what, from){
      check_num(b, what)
      if (from == 1L){
        if (length(b) != nf) stop("'", what, "' must give ", nf, " effects, one per growth factor", where, ".")
        return(matrix(b, nf, 1L))
      }
      b <- as.matrix(b)
      if (!identical(dim(b), c(nf, from))) stop("'", what, "' must be a ", nf, " x ", from, " matrix", where, ".")
      if (any(b[upper.tri(b)] != 0)) stop("In '", what, "', a growth factor affects only itself and later ones", where, ".")
      b
    }
    p$beta_xm <- path(p$beta_xm, "beta_xm", if (long) nf else 1L)
    p$beta_xy <- path(p$beta_xy, "beta_xy", if (long) nf else 1L)
    p$beta_my <- path(p$beta_my, "beta_my", nf)
    p$residual <- check_cov(p$residual, if (long) 3L else 2L, "residual", if (long) "x, m, y" else "m, y")
  }
  p
}

## n draws from a multivariate normal distribution, without packages beyond stats.
.sim_mvn <- function(n, mean, cov, what){
  cov <- as.matrix(cov)
  R <- tryCatch(chol(cov), error = function(e) stop(what, " must be positive definite.", call. = FALSE))
  sweep(matrix(rnorm(n * length(mean)), n) %*% R, 2, mean, "+")
}

## Curve without residual: each person's growth factors at their own occasions.
.sim_curve <- function(form, gf, occasions, rel_rate){
  e0 <- gf[, "eta0"]; e1 <- gf[, "eta1"]
  switch(EXPR = form,
    LIN = e0 + e1 * occasions,
    QUAD = e0 + e1 * occasions + gf[, "eta2"] * occasions^2,
    EXP = e0 + e1 * (1 - exp(-gf[, "gamma"] * occasions)),
    JB = e0 + e1 * occasions + gf[, "eta2"] * (exp(gf[, "gamma"] * occasions) - 1),
    BLS = e0 + e1 * pmin(occasions, gf[, "gamma"]) + gf[, "eta2"] * pmax(occasions - gf[, "gamma"], 0),
    NonP = e0 + e1 * .sim_cum_rates(occasions, rel_rate)
  )
}

## Cumulative change of a unit slope through intervals with relative rates.
.sim_cum_rates <- function(occasions, rel_rate){
  t(apply(occasions, 1, function(t) c(0, cumsum(rel_rate * diff(t)))))
}

## Mediation curve: linear, or a bilinear spline parameterized at the knot.
.sim_med_curve <- function(form, gf, gamma, occasions){
  if (form == "LIN") return(gf[, "eta0"] + gf[, "eta1"] * occasions)
  gf[, "etar"] + gf[, "eta1"] * pmin(0, occasions - gamma) + gf[, "eta2"] * pmax(0, occasions - gamma)
}

## Residuals of several processes at each occasion: one nk x nT matrix per process.
.sim_residuals <- function(nk, nT, cov){
  E <- .sim_mvn(nk * nT, rep(0, nrow(cov)), cov, "the residual covariance matrix")
  lapply(seq_len(ncol(E)), function(q) matrix(E[, q], nk, nT))
}

.sim_class_uni <- function(p, occasions, ctx, k){
  nk <- nrow(occasions); factors <- ctx$factors; random <- ctx$random
  gf <- matrix(p$mean[factors], nk, length(factors), byrow = TRUE, dimnames = list(NULL, factors))
  zeta <- .sim_mvn(nk, rep(0, length(random)), p$cov, paste0("the 'cov' of class ", k))
  cols <- list()
  if (!is.null(p[["beta"]])){
    x <- .sim_mvn(nk, p$tic_mean, p$tic_cov, paste0("the 'tic_cov' of class ", k))
    cols$.TIC <- x
    zeta <- zeta + sweep(x, 2, p$tic_mean) %*% t(p[["beta"]])
  }
  gf[, random] <- gf[, random] + zeta
  cols$.y1 <- .sim_curve(ctx$form, gf, occasions, p$rel_rate) +
    matrix(rnorm(nk * ctx$nT, sd = sqrt(p$residual)), nk, ctx$nT)
  list(columns = cols, gf = gf)
}

.sim_class_mgm <- function(p, occasions, ctx, k){
  nk <- nrow(occasions); factors <- ctx$factors; random <- ctx$random; nr <- length(random)
  zeta <- .sim_mvn(nk, rep(0, ncol(p$cov)), p$cov, paste0("the 'cov' of class ", k))
  E <- .sim_residuals(nk, ctx$nT, p$residual)
  cols <- list(); gfs <- list()
  for (q in seq_along(ctx$y_var)){
    P <- ctx$y_var[q]
    gf <- matrix(p$mean[[P]][factors], nk, length(factors), byrow = TRUE, dimnames = list(NULL, factors))
    gf[, random] <- gf[, random] + zeta[, (q - 1) * nr + seq_len(nr), drop = FALSE]
    cols[[paste0(".y", q)]] <- .sim_curve(ctx$form, gf, occasions, p$rel_rate[[P]]) + E[[q]]
    gfs[[q]] <- gf; colnames(gfs[[q]]) <- paste0(P, "_", factors)
  }
  list(columns = cols, gf = do.call(cbind, gfs))
}

.sim_class_tvc <- function(p, occasions, ctx, k){
  nk <- nrow(occasions); nT <- ctx$nT; factors <- ctx$factors; random <- ctx$random
  nTIC <- if (is.null(p[["beta"]])) 0L else ncol(p[["beta"]])
  gf <- matrix(p$mean[factors], nk, length(factors), byrow = TRUE, dimnames = list(NULL, factors))
  zeta <- .sim_mvn(nk, rep(0, length(random)), p$cov, paste0("the 'cov' of class ", k))
  cols <- list()
  if (ctx$decompose == 0){
    if (nTIC > 0){
      x <- .sim_mvn(nk, p$tic_mean, p$tic_cov, paste0("the 'tic_cov' of class ", k))
      cols$.TIC <- x
      zeta <- zeta + sweep(x, 2, p$tic_mean) %*% t(p[["beta"]])
    }
    gf[, random] <- gf[, random] + zeta
    tvc <- matrix(rnorm(nk * nT, mean = rep(p$tvc$mean, each = nk), sd = rep(sqrt(p$tvc$var), each = nk)), nk, nT)
    y <- .sim_curve(ctx$form, gf, occasions, p$rel_rate) + p$kappa * tvc +
      matrix(rnorm(nk * nT, sd = sqrt(p$residual)), nk, nT)
    xgf <- NULL
  } else{
    # TICs and the TVC's latent baseline and slope, jointly
    mean <- c(p$tic_mean, p$tvc$mean)
    S <- matrix(0, nTIC + 2L, nTIC + 2L)
    if (nTIC > 0){
      S[seq_len(nTIC), seq_len(nTIC)] <- p$tic_cov
      S[seq_len(nTIC), nTIC + 1L] <- S[nTIC + 1L, seq_len(nTIC)] <- p$cov_bl
    }
    S[nTIC + 1:2, nTIC + 1:2] <- p$tvc$cov
    W <- .sim_mvn(nk, mean, S, paste0("the joint covariance of the TICs and the TVC growth factors of class ", k))
    if (nTIC > 0){
      cols$.TIC <- W[, seq_len(nTIC), drop = FALSE]
      zeta <- zeta + sweep(cols$.TIC, 2, p$tic_mean) %*% t(p[["beta"]])
    }
    x0 <- W[, nTIC + 1L]; x1 <- W[, nTIC + 2L]
    zeta <- zeta + outer(x0 - p$tvc$mean[["eta0"]], p$beta_tvc)
    gf[, random] <- gf[, random] + zeta
    lag <- t(apply(occasions, 1, diff))
    slope <- outer(x1, p$tvc$rel_rate)
    true_x <- x0 + cbind(0, t(apply(slope * lag, 1, cumsum)))
    state <- switch(EXPR = as.character(ctx$decompose), "1" = cbind(0, slope), "2" = cbind(0, slope * lag),
                    "3" = true_x - x0)
    E <- .sim_residuals(nk, nT, matrix(c(p$tvc$residual, p$res_cov, p$res_cov, p$residual), 2))
    tvc <- true_x + E[[1]]
    y <- .sim_curve(ctx$form, gf, occasions, p$rel_rate) + p$kappa * state + E[[2]]
    xgf <- cbind(TVC_eta0 = x0, TVC_eta1 = x1)
  }
  cols$.y1 <- y
  cols$.TVC <- tvc
  colnames(gf) <- paste0("Y_", factors)
  list(columns = cols, gf = cbind(gf, xgf))
}

.sim_class_med <- function(p, occasions, ctx, k){
  nk <- nrow(occasions); nT <- ctx$nT; random <- ctx$random; form <- ctx$form
  long <- ctx$x_type == "longitudinal"
  draw <- function(v) .sim_mvn(nk, rep(0, length(random)), p[[v]]$cov, paste0("the '", v, "$cov' of class ", k))
  cols <- list()
  if (long){
    xg <- sweep(draw("x"), 2, p$x$mean[random], "+"); colnames(xg) <- random
    dx <- sweep(xg, 2, p$x$mean[random])
  } else{
    xb <- rnorm(nk, p$x$mean, sqrt(p$x$var))
    dx <- matrix(xb - p$x$mean, nk, 1L)
    cols$.XB <- matrix(xb, nk, 1L)
  }
  mg <- sweep(draw("m") + dx %*% t(p$beta_xm), 2, p$m$mean[random], "+"); colnames(mg) <- random
  yg <- sweep(draw("y") + dx %*% t(p$beta_xy) + sweep(mg, 2, p$m$mean[random]) %*% t(p$beta_my), 2,
              p$y$mean[random], "+")
  colnames(yg) <- random
  E <- .sim_residuals(nk, nT, p$residual)
  if (long){
    cols$.X <- .sim_med_curve(form, xg, p$x$mean[["gamma"]], occasions) + E[[1]]
    E <- E[-1]
  }
  cols$.M <- .sim_med_curve(form, mg, p$m$mean[["gamma"]], occasions) + E[[1]]
  cols$.y1 <- .sim_med_curve(form, yg, p$y$mean[["gamma"]], occasions) + E[[2]]
  colnames(mg) <- paste0("M_", random); colnames(yg) <- paste0("Y_", random)
  if (long) colnames(xg) <- paste0("X_", random)
  list(columns = cols, gf = if (long) cbind(xg, mg, yg) else cbind(mg, yg))
}

## Lower-triangle (i <= j) pairs of labels, in the order of the estimate names.
.sim_pairs <- function(labels){
  idx <- which(upper.tri(diag(length(labels)), diag = TRUE), arr.ind = TRUE)
  idx[order(idx[, "row"], idx[, "col"]), , drop = FALSE]
}

## Candidate true values of one class, by estimate name; .sim_truth keeps those that the fit reports.
.sim_candidates <- function(p, ctx){
  v <- c(); lab <- .sim_label
  cf <- ctx$curveFun; intrinsic <- ctx$intrinsic; factors <- ctx$factors; random <- ctx$random
  base <- setdiff(factors, "gamma")
  gamma_names <- function(prefix, value){
    out <- c()
    for (nm in c(.gamma_noTIC_label(cf, FALSE), .gamma_noTIC_label(cf, TRUE), .gamma_mueta_label(cf, TRUE))){
      if (!is.null(nm)) out[paste0(prefix, nm)] <- value
    }
    out
  }
  psi <- function(prefix, cov, labels = lab(random), suffix = ""){
    idx <- .sim_pairs(labels); out <- c()
    out[paste0(prefix, labels[idx[, "row"]], labels[idx[, "col"]], suffix)] <- cov[idx]
    out
  }
  if (ctx$sub_Model %in% c("LGCM", "TVC")){
    pre <- if (ctx$sub_Model == "TVC") "Y_" else ""
    nTIC <- if (is.null(p[["beta"]])) 0L else ncol(p[["beta"]])
    shift <- if (nTIC > 0) drop(p[["beta"]] %*% p$tic_mean) else rep(0, length(random))
    if (ctx$sub_Model == "TVC" && ctx$decompose != 0) shift <- shift + p$beta_tvc * p$tvc$mean[["eta0"]]
    alpha <- p$mean[random] - shift
    v[paste0(pre, "alpha", lab(base))] <- alpha[base]
    v[paste0(pre, "mueta", lab(base))] <- p$mean[base]
    if ("gamma" %in% factors){
      v <- c(v, gamma_names(pre, p$mean[["gamma"]]))
      v[paste0(pre, .gamma_alpha_label(cf, intrinsic))] <- if (intrinsic) alpha[["gamma"]] else p$mean[["gamma"]]
    }
    v <- c(v, psi(paste0(pre, "psi"), p$cov))
    if (ctx$form == "NonP") v[paste0(pre, "rel_rate", 2:(ctx$nT - 1))] <- p$rel_rate[2:(ctx$nT - 1)]
    v[paste0(pre, "residuals")] <- p$residual
    if (nTIC > 0){
      tic <- if (ctx$sub_Model == "TVC") "TIC" else ""
      ti <- .sim_pairs(seq_len(nTIC))
      if (ctx$sub_Model == "TVC" && nTIC == 1L){
        v[paste0("betaTIC", lab(random))] <- p[["beta"]][, 1]
        v["muTIC"] <- p$tic_mean; v["phiTIC"] <- p$tic_cov[1, 1]
      } else{
        for (t in seq_len(nTIC)) v[paste0("beta", tic, t, lab(random))] <- p[["beta"]][, t]
        v[paste0("mu", if (nzchar(tic)) tic else "x", seq_len(nTIC))] <- p$tic_mean
        v[paste0("phi", tic, ti[, "row"], ti[, "col"])] <- p$tic_cov[ti]
      }
    }
    if (ctx$sub_Model == "TVC"){
      v["kappa"] <- p$kappa
      if (ctx$decompose == 0){
        v[paste0("TVC_m", seq_len(ctx$nT))] <- p$tvc$mean
        v[paste0("TVC_v", seq_len(ctx$nT))] <- p$tvc$var
      } else{
        v[c("X_mueta0", "X_mueta1")] <- p$tvc$mean[c("eta0", "eta1")]
        v[c("X_psi00", "X_psi01", "X_psi11")] <- p$tvc$cov[c(1, 3, 4)]
        v[paste0("X_rel_rate", 2:(ctx$nT - 1))] <- p$tvc$rel_rate[2:(ctx$nT - 1)]
        v["X_residuals"] <- p$tvc$residual
        v[paste0("betaTVC", lab(random))] <- p$beta_tvc
        v["Cov_XYres"] <- p$res_cov
        if (nTIC == 1L) v["covBL"] <- p$cov_bl else if (nTIC > 1L) v[paste0("covBL", seq_len(nTIC))] <- p$cov_bl
      }
    }
    return(v)
  }
  if (ctx$sub_Model == "MGM"){
    P <- ctx$y_var; nr <- length(random)
    for (q in seq_along(P)){
      m <- p$mean[[P[q]]]
      v[paste0(P[q], "_mueta", lab(base))] <- m[base]
      if ("gamma" %in% factors) v <- c(v, gamma_names(paste0(P[q], "_"), m[["gamma"]]))
      i <- (q - 1) * nr + seq_len(nr)
      v <- c(v, psi(paste0(P[q], "_psi"), p$cov[i, i, drop = FALSE]))
      v[paste0(P[q], "_res")] <- p$residual[q, q]
      if (ctx$form == "NonP") v[paste0(P[q], "_rel_rate", 2:(ctx$nT - 1))] <- p$rel_rate[[P[q]]][2:(ctx$nT - 1)]
    }
    for (a in seq_along(P)) for (b in seq_along(P)) if (a < b){
      ia <- (a - 1) * nr + seq_len(nr); ib <- (b - 1) * nr + seq_len(nr)
      block <- p$cov[ia, ib, drop = FALSE]
      grid <- expand.grid(i = seq_len(nr), j = seq_len(nr))
      v[paste0(P[a], P[b], "_psi", lab(random)[grid$i], lab(random)[grid$j])] <- block[cbind(grid$i, grid$j)]
      v[paste0(P[a], P[b], "_res")] <- p$residual[a, b]
    }
    return(v)
  }
  # Mediation
  lr <- lab(random); long <- ctx$x_type == "longitudinal"
  mux <- if (long) p$x$mean[random] else p$x$mean
  alpha_m <- p$m$mean[random] - drop(p$beta_xm %*% mux)
  alpha_y <- p$y$mean[random] - drop(p$beta_xy %*% mux) - drop(p$beta_my %*% p$m$mean[random])
  paths <- function(b, from, to){
    out <- c()
    for (a in seq_along(lr)) for (j in seq_along(lr)) if (a <= j) out[paste0("beta", from, lr[a], to, lr[j])] <- b[j, a]
    out
  }
  if (long){
    v[paste0("muetaX", lr)] <- p$x$mean[random]
    if ("gamma" %in% factors) v["mugX"] <- p$x$mean[["gamma"]]
    v <- c(v, psi("psiX", p$x$cov, suffix = ""))
    v <- c(v, paths(p$beta_xm, "X", "M"), paths(p$beta_xy, "X", "Y"))
  } else{
    v["muX"] <- p$x$mean; v["phi11"] <- p$x$var
    v[paste0("betaM", lr)] <- p$beta_xm[, 1]
    v[paste0("betaY", lr)] <- p$beta_xy[, 1]
  }
  for (w in c("M", "Y")){
    pp <- p[[tolower(w)]]
    v[paste0("alpha", w, lr)] <- if (w == "M") alpha_m else alpha_y
    if ("gamma" %in% factors) v[paste0("mug", w)] <- pp$mean[["gamma"]]
    v[paste0("mueta", w, lr)] <- pp$mean[random]
    idx <- .sim_pairs(lr)
    v[paste0("psi", w, lr[idx[, "row"]], w, lr[idx[, "col"]], "_r")] <- pp$cov[idx]
  }
  if (long){
    idx <- .sim_pairs(lr)
    v <- v[!names(v) %in% paste0("psiX", lr[idx[, "row"]], lr[idx[, "col"]])]
    v[paste0("psiX", lr[idx[, "row"]], "X", lr[idx[, "col"]])] <- p$x$cov[idx]
  }
  v <- c(v, paths(p$beta_my, "M", "Y"))
  R <- p$residual
  if (long){
    v[c("residualsX", "residualsM", "residualsY", "residualsMX", "residualsYX", "residualsYM")] <-
      c(R[1, 1], R[2, 2], R[3, 3], R[2, 1], R[3, 1], R[3, 2])
  } else{
    v[c("residualsM", "residualsY", "residualsYM")] <- c(R[1, 1], R[2, 2], R[2, 1])
  }
  v
}

## Estimate names of one class, without the derived quantities.
.sim_fit_names <- function(ctx){
  nT <- ctx$nT
  out <- switch(EXPR = ctx$sub_Model,
    LGCM = if (ctx$y_model == "LCSM") .auto_names_LCSM(ctx$curveFun, ctx$intrinsic, ctx$growth_TIC, nT) else
      .auto_names_LGCM(ctx$curveFun, ctx$intrinsic, ctx$growth_TIC),
    TVC = .auto_names_TVC(ctx$curveFun, ctx$intrinsic, seq_len(nT), ctx$y_model, ctx$decompose, ctx$growth_TIC),
    MGM = .auto_names_MGM(ctx$curveFun, ctx$intrinsic, ctx$y_var, rep(list(seq_len(nT)), length(ctx$y_var)), ctx$y_model),
    MED = .auto_names_MED(ctx$curveFun, ctx$x_type, "X", "M", "Y"))
  out[!grepl("_est[0-9]*$|^X_abs_rate|^[Mm]ediator|^total", out)]
}

## True values named as in the Estimates of the matching fit.
.sim_truth <- function(ctx, params, prop, cluster){
  keep <- .sim_fit_names(ctx)
  one <- function(p){
    v <- .sim_candidates(p, ctx)
    missing <- setdiff(keep, names(v))
    if (length(missing)) stop("No true value for ", paste(missing, collapse = ", "), call. = FALSE)
    v[keep]
  }
  if (length(params) == 1L) return(one(params[[1]]))
  out <- unlist(lapply(seq_along(params), function(k){ v <- one(params[[k]]); stats::setNames(v, paste0("c", k, names(v))) }))
  if (!is.null(cluster)){
    beta <- cluster$beta[-1, , drop = FALSE]
    out[.cluster_tic_param_names(length(params), seq_len(ncol(beta) - 1L))] <- c(t(beta))
  } else{
    prop <- if (is.null(prop)) rep(1 / length(params), length(params)) else prop / sum(prop)
    out[paste0("p", seq_along(params)[-1])] <- log(prop[-1] / prop[1])
  }
  out
}
