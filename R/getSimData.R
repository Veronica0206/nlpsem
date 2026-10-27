#' @title Simulate Longitudinal Data from a Latent Growth Curve or Latent Change Score Model
#'
#' @description This function generates wide-format longitudinal data from a latent growth curve model or a
#' latent change score model with one of the functional forms that \code{getLGCM()} and \code{getLCSM()} fit,
#' optionally with time-invariant covariates (TICs) and several latent classes. The true parameter values are
#' attached to the data under the names that the matching fitting function uses in its \code{Estimates}, so
#' that a fit can be compared with the values that generated the data.
#'
#' @param n An integer: the number of persons. Default is \code{500}.
#' @param sub_Model A string: \code{"LGCM"} (latent growth curve model) or \code{"LCSM"} (latent change score model).
#' Default is \code{"LGCM"}.
#' @param curveFun A string specifying the functional form. Supported options for \code{sub_Model = "LGCM"} are
#' \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or
#' \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}); for
#' \code{sub_Model = "LCSM"}, \code{"quadratic"}, \code{"negative exponential"}, \code{"Jenss-Bayley"}, and
#' \code{"nonparametric"} (or \code{"NonP"}).
#' @param intrinsic A logical flag: whether the growth ratio (negative exponential and Jenss-Bayley) or the knot
#' (bilinear spline) varies between persons. Default is \code{FALSE}.
#' @param nT An integer: the number of measurement occasions. Default is \code{10}.
#' @param times A numeric vector of the \code{nT} planned occasions. Default is \code{0, 1, ..., nT - 1}.
#' @param jitter A non-negative number: each person's interior occasions are drawn uniformly within
#' \code{jitter} of the planned ones (individually varying measurement occasions); the first and the last
#' occasions are the planned ones. Default is \code{0.25}.
#' @param params A list of parameters for one class, or a list of such lists, one per latent class. The elements
#' of a class list are:
#' \itemize{
#'   \item \code{mean}: a named numeric vector of the means of the growth factors: \code{eta0} (intercept),
#'   \code{eta1} (slope; for the negative exponential, the change to the asymptote; for the bilinear spline, the
#'   slope before the knot; for the nonparametric form, the slope of the first interval), \code{eta2}
#'   (quadratic coefficient; for the Jenss-Bayley, the coefficient of the exponential term; for the bilinear
#'   spline, the slope after the knot), and \code{gamma} (growth ratio or knot), as the form requires;
#'   \item \code{cov}: the covariance matrix of the growth factors that vary between persons (in the order
#'   above, \code{gamma} last), not explained by the growth TICs;
#'   \item \code{residual}: the residual variance;
#'   \item \code{rel_rate}: for the nonparametric form, the relative rates of the \code{nT - 1} intervals; the
#'   first is \code{1};
#'   \item \code{tic_mean}, \code{tic_cov}, \code{beta}: for growth TICs, their means and covariance matrix, and
#'   the matrix of their effects on the growth factors that vary between persons (one row per growth factor,
#'   one column per TIC).
#' }
#' Default is \code{NULL}: one class with the settings of the package's simulation studies.
#' @param prop A numeric vector of class proportions for several classes without cluster TICs. Default is
#' \code{NULL} (equal proportions).
#' @param cluster A list for several classes whose membership depends on cluster TICs: \code{mean} and
#' \code{cov} of the cluster TICs, and \code{beta}, a matrix of multinomial-logit coefficients with one row per
#' class (the first row zeros, for the reference class) and one column for the intercept and each TIC.
#' Default is \code{NULL}.
#' @param growth_TIC A character vector of column names for the growth TICs. Default is \code{"ex1"},
#' \code{"ex2"}, ... when \code{params} has growth-TIC effects.
#' @param cluster_TIC A character vector of column names for the cluster TICs. Default is \code{"gx1"},
#' \code{"gx2"}, ... when \code{cluster} is given.
#' @param y_var,t_var Strings: the prefixes of the outcome and occasion columns. Defaults are \code{"Y"} and
#' \code{"T"}.
#' @param seed An integer seed for reproducible data, or \code{NULL}. The random number generator state of the
#' session is restored afterwards. Default is \code{NULL}.
#'
#' @return A data frame with columns \code{id}, \code{y_var1}, ..., \code{t_var1}, ..., the TICs and, for several
#' classes, \code{class}. Its attribute \code{"truth"} holds the true parameter values, named as in the
#' \code{Estimates} of \code{getLGCM()}/\code{getLCSM()} (one class) or \code{getMIX()} (several classes, with the
#' class weights as log odds \code{p2}, ... against the first class, or as cluster-TIC coefficients
#' \code{beta20}, ...); attribute \code{"growth_factors"} holds each person's growth factors, and attribute
#' \code{"settings"} the arguments.
#'
#' @details Growth factors and growth TICs are drawn from multivariate normal distributions, each person's
#' outcome follows the functional form at their own occasions, and normal residuals are added. With
#' \code{sub_Model = "LCSM"}, a parametric form generates the same trajectories as with \code{"LGCM"} (a latent
#' change score model describes such a trajectory by its rates of change), and the nonparametric form
#' generates change scores: the slope of each interval is \code{eta1} times the interval's relative rate.
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
#'
#' @importFrom stats rnorm runif
#' @export
getSimData <- function(n = 500, sub_Model = "LGCM", curveFun, intrinsic = FALSE, nT = 10, times = NULL,
                       jitter = 0.25, params = NULL, prop = NULL, cluster = NULL, growth_TIC = NULL,
                       cluster_TIC = NULL, y_var = "Y", t_var = "T", seed = NULL){
  validate_curveFun(curveFun)
  if (!isTRUE(sub_Model %in% c("LGCM", "LCSM"))){
    stop("sub_Model must be 'LGCM' or 'LCSM'.")
  }
  validate_curveFun_ymodel(curveFun, sub_Model)
  validate_intrinsic(intrinsic, curveFun, model_type = sub_Model)
  form <- .sim_form(curveFun)
  if (sub_Model == "LGCM" && form == "NonP"){
    stop("The nonparametric functional form is a latent change score model; use sub_Model = 'LCSM'.")
  }
  if (!is.numeric(n) || length(n) != 1L || n < 2 || n != round(n)) stop("'n' must be an integer of at least 2.")
  if (!is.numeric(nT) || length(nT) != 1L || nT < 3 || nT != round(nT)) stop("'nT' must be an integer of at least 3.")
  if (is.null(times)) times <- seq_len(nT) - 1
  if (!is.numeric(times) || length(times) != nT || any(diff(times) <= 0)){
    stop("'times' must be ", nT, " increasing numbers.")
  }
  if (!is.numeric(jitter) || length(jitter) != 1L || jitter < 0 || (nT > 2 && jitter >= min(diff(times)) / 2)){
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
  factors <- .sim_factors(form)
  random <- if (intrinsic || !"gamma" %in% factors) factors else setdiff(factors, "gamma")
  if (is.null(params)) params <- list(.sim_default_params(form, intrinsic, nT))
  if (!is.null(params$mean)) params <- list(params)
  nClass <- length(params)
  params <- lapply(seq_len(nClass), function(k) .sim_check_params(params[[k]], k, form, factors, random, nT))
  nTIC <- unique(vapply(params, function(p) if (is.null(p$beta)) 0L else ncol(p$beta), integer(1)))
  if (length(nTIC) != 1L) stop("All classes need the same number of growth TICs.")
  if (nTIC > 0 && is.null(growth_TIC)) growth_TIC <- paste0("ex", seq_len(nTIC))
  if (length(growth_TIC) != nTIC) stop("'growth_TIC' must name the ", nTIC, " growth TIC(s) in 'params'.")

  # Class membership
  X <- NULL
  if (nClass == 1L){
    if (!is.null(prop) || !is.null(cluster)) stop("'prop' and 'cluster' need parameters for several classes.")
    class <- rep(1L, n)
  } else if (!is.null(cluster)){
    if (!is.null(prop)) stop("Give either 'prop' or 'cluster', not both.")
    nX <- length(cluster$mean)
    if (!is.matrix(cluster$beta) || !identical(dim(cluster$beta), c(nClass, nX + 1L)) || any(cluster$beta[1, ] != 0)){
      stop("'cluster$beta' must be a ", nClass, " x ", nX + 1L, " matrix with a first row of zeros.")
    }
    if (is.null(cluster_TIC)) cluster_TIC <- paste0("gx", seq_len(nX))
    if (length(cluster_TIC) != nX) stop("'cluster_TIC' must name the ", nX, " cluster TIC(s).")
    X <- .sim_mvn(n, cluster$mean, cluster$cov, "cluster$cov")
    eta <- cbind(1, X) %*% t(cluster$beta)
    prob <- exp(eta - apply(eta, 1, max)); prob <- prob / rowSums(prob)
    class <- 1L + rowSums(runif(n) > t(apply(prob, 1, cumsum))[, -nClass, drop = FALSE])
  } else{
    if (is.null(prop)) prop <- rep(1 / nClass, nClass)
    if (!is.numeric(prop) || length(prop) != nClass || any(prop <= 0)) stop("'prop' must be ", nClass, " positive numbers.")
    prop <- prop / sum(prop)
    class <- sample.int(nClass, n, replace = TRUE, prob = prop)
  }

  # Individual occasions
  occasions <- matrix(times, n, nT, byrow = TRUE)
  if (jitter > 0) occasions[, 2:(nT - 1)] <- occasions[, 2:(nT - 1)] + runif(n * (nT - 2), -jitter, jitter)

  # Growth factors, growth TICs and outcomes, class by class
  Y <- matrix(NA_real_, n, nT)
  TIC <- if (nTIC > 0) matrix(NA_real_, n, nTIC) else NULL
  GF <- matrix(NA_real_, n, length(factors), dimnames = list(NULL, factors))
  for (k in seq_len(nClass)){
    rows <- which(class == k)
    if (!length(rows)) next
    p <- params[[k]]
    gf <- matrix(p$mean[factors], length(rows), length(factors), byrow = TRUE, dimnames = list(NULL, factors))
    zeta <- .sim_mvn(length(rows), rep(0, length(random)), p$cov, paste0("the 'cov' of class ", k))
    if (nTIC > 0){
      x <- .sim_mvn(length(rows), p$tic_mean, p$tic_cov, paste0("the 'tic_cov' of class ", k))
      TIC[rows, ] <- x
      zeta <- zeta + sweep(x, 2, p$tic_mean) %*% t(p$beta)
    }
    gf[, random] <- gf[, random] + zeta
    GF[rows, ] <- gf
    Y[rows, ] <- .sim_curve(form, gf, occasions[rows, , drop = FALSE], p$rel_rate) +
      matrix(rnorm(length(rows) * nT, sd = sqrt(p$residual)), length(rows), nT)
  }

  dat <- data.frame(id = seq_len(n))
  for (j in seq_len(nT)) dat[[paste0(y_var, j)]] <- Y[, j]
  for (j in seq_len(nT)) dat[[paste0(t_var, j)]] <- occasions[, j]
  if (nTIC > 0) for (j in seq_len(nTIC)) dat[[growth_TIC[j]]] <- TIC[, j]
  if (!is.null(X)) for (j in seq_along(cluster_TIC)) dat[[cluster_TIC[j]]] <- X[, j]
  if (nClass > 1L) dat$class <- class

  attr(dat, "truth") <- .sim_truth(sub_Model, curveFun, intrinsic, form, factors, random, params, growth_TIC,
                                   nT, prop, cluster)
  attr(dat, "growth_factors") <- GF
  attr(dat, "settings") <- list(n = n, sub_Model = sub_Model, curveFun = curveFun, intrinsic = intrinsic, nT = nT,
                                times = times, jitter = jitter, params = params, prop = prop, cluster = cluster,
                                growth_TIC = growth_TIC, cluster_TIC = cluster_TIC, y_var = y_var, t_var = t_var,
                                seed = seed)
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

## Settings of the package's simulation studies: intercept variance 25, other variances 1 (the growth ratio
## 0.1^2, the knot 0.3^2), correlations 0.3 and residual variance 1.
.sim_default_params <- function(form, intrinsic, nT){
  d <- list(LIN = list(mean = c(eta0 = 50, eta1 = 5), sd = c(5, 1)),
            QUAD = list(mean = c(eta0 = 50, eta1 = 20, eta2 = -1), sd = c(5, 1, 0.3)),
            EXP = list(mean = c(eta0 = 50, eta1 = 30, gamma = 0.6), sd = c(5, 3, 0.1)),
            JB = list(mean = c(eta0 = 50, eta1 = 2.5, eta2 = -30, gamma = -0.7), sd = c(5, 1, 3, 0.1)),
            BLS = list(mean = c(eta0 = 100, eta1 = -5, eta2 = -1.8, gamma = 4.5), sd = c(5, 1, 1, 0.3)),
            NonP = list(mean = c(eta0 = 50, eta1 = 5), sd = c(5, 1)))[[form]]
  sd <- d$sd
  if (!intrinsic && form %in% c("EXP", "JB", "BLS")) sd <- sd[-length(sd)]
  R <- matrix(0.3, length(sd), length(sd)); diag(R) <- 1
  out <- list(mean = d$mean, cov = R * outer(sd, sd), residual = 1)
  if (form == "NonP") out$rel_rate <- seq(1, by = -0.1 * 9 / (nT - 1), length.out = nT - 1)
  out
}

.sim_check_params <- function(p, k, form, factors, random, nT){
  where <- paste0(" (class ", k, ")")
  if (!is.list(p) || is.null(p$mean) || is.null(p$cov) || is.null(p$residual)){
    stop("Each class needs 'mean', 'cov' and 'residual'", where, ".")
  }
  if (is.null(names(p$mean)) || !all(factors %in% names(p$mean))){
    stop("'mean' must name ", paste(factors, collapse = ", "), where, ".")
  }
  p$cov <- as.matrix(p$cov)
  if (!identical(dim(p$cov), rep(length(random), 2L))){
    stop("'cov' must be a ", length(random), " x ", length(random), " matrix for ", paste(random, collapse = ", "), where, ".")
  }
  if (!is.numeric(p$residual) || length(p$residual) != 1L || p$residual <= 0) stop("'residual' must be positive", where, ".")
  if (form == "NonP"){
    if (is.null(p$rel_rate)) stop("The nonparametric form needs 'rel_rate'", where, ".")
    if (length(p$rel_rate) != nT - 1L || p$rel_rate[1] != 1){
      stop("'rel_rate' must give the ", nT - 1L, " intervals, the first 1", where, ".")
    }
  }
  if (!is.null(p$beta)){
    p$beta <- as.matrix(p$beta)
    nTIC <- ncol(p$beta)
    if (nrow(p$beta) != length(random) || is.null(p$tic_mean) || length(p$tic_mean) != nTIC ||
        !identical(dim(as.matrix(p$tic_cov)), c(nTIC, nTIC))){
      stop("Growth TICs need 'tic_mean', 'tic_cov' and a 'beta' with one row per growth factor that varies",
           " between persons", where, ".")
    }
    p$tic_cov <- as.matrix(p$tic_cov)
  }
  p
}

## n draws from a multivariate normal distribution, without packages beyond stats.
.sim_mvn <- function(n, mean, cov, what){
  cov <- as.matrix(cov)
  R <- tryCatch(chol(cov), error = function(e) stop(what, " must be positive definite.", call. = FALSE))
  sweep(matrix(rnorm(n * length(mean)), n) %*% R, 2, mean, "+")
}

## Outcome without residual: each person's curve at their own occasions.
.sim_curve <- function(form, gf, occasions, rel_rate){
  e0 <- gf[, "eta0"]; e1 <- gf[, "eta1"]
  switch(EXPR = form,
    LIN = e0 + e1 * occasions,
    QUAD = e0 + e1 * occasions + gf[, "eta2"] * occasions^2,
    EXP = e0 + e1 * (1 - exp(-gf[, "gamma"] * occasions)),
    JB = e0 + e1 * occasions + gf[, "eta2"] * (exp(gf[, "gamma"] * occasions) - 1),
    BLS = e0 + e1 * pmin(occasions, gf[, "gamma"]) + gf[, "eta2"] * pmax(occasions - gf[, "gamma"], 0),
    NonP = e0 + e1 * t(apply(occasions, 1, function(t) c(0, cumsum(rel_rate * diff(t)))))
  )
}

## True values named as in the Estimates of the matching fit.
.sim_truth <- function(sub_Model, curveFun, intrinsic, form, factors, random, params, growth_TIC, nT, prop, cluster){
  label <- c(eta0 = "0", eta1 = "1", eta2 = "2", gamma = "g")
  one <- function(p){
    v <- c()
    base <- setdiff(factors, "gamma")
    nTIC <- length(growth_TIC)
    gamma_name <- .gamma_noTIC_label(curveFun, intrinsic)
    if (nTIC > 0){
      alpha <- p$mean[random] - drop(p$beta %*% p$tic_mean)
      v[paste0("alpha", label[base])] <- alpha[base]
      if ("gamma" %in% factors){
        v[.gamma_alpha_label(curveFun, intrinsic)] <- if (intrinsic) alpha[["gamma"]] else p$mean[["gamma"]]
      }
    } else{
      v[paste0("mueta", label[base])] <- p$mean[base]
      if (!is.null(gamma_name)) v[gamma_name] <- p$mean[["gamma"]]
    }
    idx <- which(lower.tri(p$cov, diag = TRUE), arr.ind = TRUE)
    idx <- idx[order(idx[, "col"], idx[, "row"]), , drop = FALSE]
    v[paste0("psi", label[random[idx[, "col"]]], label[random[idx[, "row"]]])] <- p$cov[idx]
    if (form == "NonP") v[paste0("rel_rate", 2:(nT - 1))] <- p$rel_rate[2:(nT - 1)]
    v["residuals"] <- p$residual
    if (nTIC > 0){
      for (t in seq_len(nTIC)) v[paste0("beta", t, label[random])] <- p$beta[, t]
      v[paste0("mux", seq_len(nTIC))] <- p$tic_mean
      ti <- which(lower.tri(p$tic_cov, diag = TRUE), arr.ind = TRUE)
      ti <- ti[order(ti[, "col"], ti[, "row"]), , drop = FALSE]
      v[paste0("phi", ti[, "col"], ti[, "row"])] <- p$tic_cov[ti]
      v[paste0("mueta", label[base])] <- p$mean[base]
      gm <- .gamma_mueta_label(curveFun, intrinsic)
      if (!is.null(gm)) v[gm] <- p$mean[["gamma"]]
    }
    names_fit <- if (sub_Model == "LGCM") .auto_names_LGCM(curveFun, intrinsic, growth_TIC) else
      .auto_names_LCSM(curveFun, intrinsic, growth_TIC, nT)
    v[intersect(names_fit, names(v))]
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
