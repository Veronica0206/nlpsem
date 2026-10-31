#' Validate that parameter names are provided when paramOut is TRUE
#' @param paramOut Logical flag for parameter output.
#' @param names Character vector of parameter names.
#' @keywords internal
#' @noRd
validate_paramOut <- function(paramOut, names) {
  if (!is.logical(paramOut) || length(paramOut) != 1L || is.na(paramOut)) {
    stop("paramOut must be a single logical value (TRUE or FALSE).")
  }
  if (!is.null(names)) .validate_output_names(names)
  invisible(NULL)
}

# Parameter labels are a one-to-one mapping, never a recyclable data-frame
# column. The optional derived names preserve standalone LCSM's documented
# core-only naming convention; arbitrary shorter vectors are rejected.
#' @keywords internal
#' @noRd
.validate_output_names <- function(names, expected_length = NULL, derived_names = character()) {
  if (!is.character(names) || !is.null(dim(names))) {
    stop("names must be a character vector or NULL.", call. = FALSE)
  }
  if (!length(names) || anyNA(names) || any(!nzchar(trimws(names))) || anyDuplicated(names)) {
    stop("names must contain non-missing, non-empty, unique parameter names.", call. = FALSE)
  }
  if (!is.null(expected_length)) {
    if (length(derived_names) && length(names) == expected_length - length(derived_names)) {
      names <- c(names, derived_names)
    }
    if (length(names) != expected_length) {
      stop("names must contain exactly ", expected_length, " parameter names",
           if (length(derived_names)) paste0(" (or ", expected_length - length(derived_names),
                                             " core names for this LCSM)"),
           "; received ", length(names), ".", call. = FALSE)
    }
    if (anyDuplicated(names)) {
      stop("names must remain unique after derived parameter names are appended.", call. = FALSE)
    }
  }
  names
}

#' Validate res_scale values are between 0 and 1 (exclusive)
#'
#' Handles scalar, vector, and list inputs for res_scale.
#' @param res_scale A numeric value, vector, or list of scaling factors.
#' @param is_list Logical indicating whether res_scale is expected as a list
#'   (e.g., for getMGroup or getMIX). Default is \code{FALSE}.
#' @keywords internal
#' @noRd
validate_res_scale <- function(res_scale, is_list = FALSE) {
  if (is.null(res_scale)) return(invisible(NULL))
  if (is_list) {
    if (any(unlist(res_scale) <= 0 | unlist(res_scale) >= 1)) {
      stop("Please enter a value between 0 and 1 (exclusive) for res_scale!")
    }
  } else {
    if (any(res_scale <= 0) || any(res_scale >= 1)) {
      stop("Please enter a value between 0 and 1 (exclusive) for res_scale!")
    }
  }
}

#' Validate intrinsic and curveFun compatibility
#'
#' Checks that intrinsic nonlinear models are not paired with linear,
#' quadratic, or nonparametric functional forms.
#' @param intrinsic Logical flag for intrinsic nonlinear model.
#' @param curveFun A string specifying the functional form.
#' @param model_type A string identifying the model context, used to
#'   customize the error message. Supported values: \code{"LGCM"},
#'   \code{"LCSM"}, \code{"general"}, \code{"sub_Model"}.
#' @keywords internal
#' @noRd
validate_intrinsic <- function(intrinsic, curveFun, model_type = "general") {
  if (is.null(intrinsic)) return(invisible(NULL))

  if (model_type == "LGCM") {
    if (intrinsic && curveFun %in% c("linear", "LIN", "quadratic", "QUAD")) {
      stop("An intrinsic nonlinear function should be negative exponential, Jenss-Bayley, or bilinear spline for a LGCM!")
    }
  } else if (model_type == "LCSM") {
    if (intrinsic && curveFun %in% c("nonparametric", "NonP", "quadratic", "QUAD")) {
      stop("An intrinsic nonlinear function should be negative exponential or Jenss-Bayley for a LCSM!")
    }
  } else if (model_type == "sub_Model") {
    if (intrinsic && curveFun %in% c("linear", "LIN", "quadratic", "QUAD")) {
      stop("An intrinsic nonlinear function should be one of the following: negative exponential, Jenss-Bayley, or bilinear spline!")
    }
  } else {
    # general: used by getTVCmodel and getMGM
    if (intrinsic && curveFun %in% c("linear", "LIN", "quadratic", "QUAD", "nonparametric", "NonP")) {
      stop("An intrinsic nonlinear function should be negative exponential, Jenss-Bayley, or bilinear spline!")
    }
  }
}

#' Validate sub_Model-specific arguments for getMGroup and getMIX
#'
#' Checks that required arguments are provided depending on the sub_Model.
#' @param sub_Model A string specifying the sub-model.
#' @param intrinsic Logical flag for intrinsic nonlinear model.
#' @param curveFun A string specifying the functional form.
#' @param res_cor Residual correlation value(s).
#' @param y_model A string specifying the longitudinal outcome model.
#' @param x_type A string for predictor type in mediation models.
#' @param x_var Predictor variable specification.
#' @param m_var Mediator variable specification.
#' @param y_var Outcome variable specification.
#' @param decompose Decomposition option for TVC models.
#' @keywords internal
#' @noRd
validate_sub_Model <- function(sub_Model, intrinsic = NULL, curveFun = NULL,
                               res_cor = NULL, y_model = NULL, x_type = NULL,
                               x_var = NULL, m_var = NULL, y_var = NULL,
                               decompose = NULL) {
  if (sub_Model %in% c("LGCM", "LCSM", "TVC", "MGM")) {
    if (is.null(intrinsic)) {
      stop("Please specify if an intrinsic nonlinear model will be built if the sub_Model is LGCM, LCSM, MGM, or TVC!")
    } else {
      validate_intrinsic(intrinsic, curveFun, model_type = "sub_Model")
    }
  }
  # res_cor is now optional — when NULL, data-driven estimation is used in the initial value helpers
  if (sub_Model %in% c("TVC", "MGM")) {
    if (is.null(y_model)) {
      stop("Please specify y_model (either LGCM or LCSM) if the sub_Model is TVC or MGM!")
    }
  }
  if (sub_Model %in% c("MED")) {
    if (is.null(x_type)) {
      stop("Please specify x_type for a longitudinal mediation model!")
    } else if (is.null(x_var) || is.null(m_var) || is.null(y_var)) {
      stop("Please specify predictor, mediator and outcome for a longitudinal mediation model!")
    } else if (!(curveFun %in% c("linear", "LIN", "bilinear spline", "BLS"))) {
      stop("Longitudinal mediation model only allows for linear or bilinear spline functional form!")
    }
  }
}

#' Validate that curveFun is a supported value
#'
#' @param curveFun A string specifying the functional form.
#' @keywords internal
#' @noRd
validate_curveFun <- function(curveFun) {
  supported <- c("linear", "LIN", "quadratic", "QUAD", "negative exponential", "EXP",
                  "Jenss-Bayley", "JB", "bilinear spline", "BLS", "nonparametric", "NonP")
  if (!curveFun %in% supported) {
    stop("Unsupported curveFun '", curveFun, "'. Supported options: ",
         paste(supported, collapse = ", "), ".")
  }
}

#' Validate curveFun compatibility with y_model
#'
#' LCSM supports only NonP, QUAD, EXP, JB (not LIN or BLS).
#' LGCM supports LIN, QUAD, EXP, JB, BLS.
#' @param curveFun A string specifying the functional form.
#' @param y_model A string specifying the outcome model type ("LGCM" or "LCSM").
#' @keywords internal
#' @noRd
validate_curveFun_ymodel <- function(curveFun, y_model) {
  if (is.null(y_model)) return(invisible(NULL))
  if (!y_model %in% c("LGCM", "LCSM")) {
    stop("y_model must be 'LGCM' or 'LCSM', got '", y_model, "'.")
  }
  if (y_model == "LCSM" && curveFun %in% c("linear", "LIN", "bilinear spline", "BLS")) {
    stop("LCSM does not support '", curveFun, "' functional form. ",
         "Supported options for LCSM: nonparametric (NonP), quadratic (QUAD), ",
         "negative exponential (EXP), Jenss-Bayley (JB).")
  }
}

#' Validate decompose value
#'
#' @param decompose An integer specifying the decomposition option.
#' @keywords internal
#' @noRd
validate_decompose <- function(decompose) {
  if (is.null(decompose)) return(invisible(NULL))
  if (!decompose %in% c(0, 1, 2, 3)) {
    stop("decompose must be 0, 1, 2, or 3, got '", decompose, "'.")
  }
}

#' Validate x_type value
#'
#' @param x_type A string indicating the predictor type.
#' @keywords internal
#' @noRd
validate_x_type <- function(x_type) {
  if (is.null(x_type)) return(invisible(NULL))
  if (!x_type %in% c("baseline", "longitudinal")) {
    stop("x_type must be 'baseline' or 'longitudinal', got '", x_type, "'.")
  }
}

#' Validate res_cor values are between -1 and 1 (exclusive)
#'
#' Handles scalar, vector, and list inputs for res_cor.
#' @param res_cor A numeric value, vector, or list of correlation values.
#' @keywords internal
#' @noRd
validate_res_cor <- function(res_cor) {
  if (is.null(res_cor)) return(invisible(NULL))
  vals <- unlist(res_cor)
  if (any(vals <= -1 | vals >= 1)) {
    stop("Please enter a value between -1 and 1 (exclusive) for res_cor!")
  }
}

#' Validate that expected columns exist in data
#'
#' Checks whether columns constructed from variable prefixes and record indices
#' actually exist in the provided data frame. Stops with a clear error message
#' listing all missing columns if any are not found.
#'
#' @param dat A data frame.
#' @param t_var A string or character vector of time variable prefix(es).
#' @param y_var A string or character vector of outcome variable prefix(es).
#' @param records A numeric vector or list of numeric vectors specifying wave indices.
#' @param growth_TIC Optional character vector of TIC column names.
#' @param TVC Optional string specifying the TVC variable prefix.
#' @param m_var Optional string specifying the mediator variable prefix.
#' @param x_var Optional string specifying the predictor variable prefix.
#' @param x_type Optional string indicating the predictor type ("baseline" or "longitudinal").
#' @param grp_var Optional string specifying the grouping variable column name.
#' @param cluster_TIC Optional character vector of cluster TIC column names.
#' @keywords internal
#' @noRd
validate_columns <- function(dat, t_var, y_var, records,
                             growth_TIC = NULL, TVC = NULL,
                             m_var = NULL, x_var = NULL, x_type = NULL,
                             grp_var = NULL, cluster_TIC = NULL) {
  expected <- character(0)

  # Helper: generate prefix+index columns for a single prefix and record vector
  .cols_from <- function(prefix, recs) {
    paste0(prefix, recs)
  }

  # Handle records as list (for MGM/MED) or vector
  if (is.list(records)) {
    # t_var and y_var may be vectors of same length as records
    for (i in seq_along(records)) {
      tv <- if (length(t_var) >= i) t_var[i] else t_var[1]
      expected <- c(expected, .cols_from(tv, records[[i]]))
    }
    # y_var may be single or vector
    if (length(y_var) == length(records)) {
      for (i in seq_along(records)) {
        expected <- c(expected, .cols_from(y_var[i], records[[i]]))
      }
    } else {
      expected <- c(expected, .cols_from(y_var[1], records[[1]]))
    }
    # m_var uses second record set
    if (!is.null(m_var) && length(records) >= 2) {
      expected <- c(expected, .cols_from(m_var, records[[2]]))
    }
    # x_var uses last record set for longitudinal, or is just a column name for baseline
    if (!is.null(x_var)) {
      if (!is.null(x_type) && x_type == "baseline") {
        expected <- c(expected, x_var)
      } else if (length(records) >= 3) {
        expected <- c(expected, .cols_from(x_var, records[[3]]))
      }
    }
  } else {
    # Single records vector
    for (tv in t_var) {
      expected <- c(expected, .cols_from(tv, records))
    }
    for (yv in y_var) {
      expected <- c(expected, .cols_from(yv, records))
    }
    if (!is.null(m_var)) {
      expected <- c(expected, .cols_from(m_var, records))
    }
    if (!is.null(x_var)) {
      if (!is.null(x_type) && x_type == "baseline") {
        expected <- c(expected, x_var)
      } else {
        expected <- c(expected, .cols_from(x_var, records))
      }
    }
    if (!is.null(TVC)) {
      expected <- c(expected, .cols_from(TVC, records))
    }
  }

  # TICs and grouping vars
  if (!is.null(growth_TIC)) expected <- c(expected, growth_TIC)
  if (!is.null(grp_var)) expected <- c(expected, grp_var)
  if (!is.null(cluster_TIC)) expected <- c(expected, cluster_TIC)

  # De-duplicate and check
  expected <- unique(expected)
  missing_cols <- setdiff(expected, colnames(dat))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in data: ", paste(missing_cols, collapse = ", "))
  }
  invisible(NULL)
}

#' Recode a grouping variable to consecutive integer codes 1..K
#'
#' Multiple-group helpers select each group with \code{dat[, grp_var] == k}, so the codes must be
#' \code{1, ..., K}. Distinct non-missing values are sorted and mapped to their rank.
#' @param x The grouping variable (numeric, character, or factor).
#' @param grp_var Column name, used in messages.
#' @return A list with \code{codes} (integer, \code{NA} where \code{x} is missing), \code{nClass},
#'   and \code{levels} (the original values in the order of the new codes).
#' @keywords internal
#' @noRd
.recode_group_codes <- function(x, grp_var = "grp_var") {
  missing <- is.na(x)
  levels <- sort(unique(x[!missing]))
  if (length(levels) < 2) {
    stop("The grouping variable '", grp_var, "' must have at least two distinct non-missing values.")
  }
  codes <- match(x, levels)
  if (any(missing)) {
    warning(sum(missing), " row(s) with a missing value in '", grp_var,
            "' were removed before fitting the multiple-group model.", call. = FALSE)
  }
  already_consecutive <- is.numeric(levels) &&
    isTRUE(all.equal(as.numeric(levels), as.numeric(seq_along(levels))))
  if (!already_consecutive) {
    message("Groups in '", grp_var, "' are recoded for estimation: ",
            paste(sprintf("%s -> %d", as.character(levels), seq_along(levels)), collapse = ", "))
  }
  list(codes = as.integer(codes), nClass = length(levels), levels = levels)
}

#' Check that the process prefixes of a multivariate model give unambiguous parameter labels
#'
#' Between-process parameters are labelled by joining the prefixes of the two processes (the pair
#' \code{"A"}, \code{"B"} gives \code{"AB_psi00"}), so a joined pair must not equal a process prefix or
#' another pair's; latent change scores use the prefixes in lower case, so these must differ too.
#'
#' @param y_var Character vector of process prefixes.
#' @param y_model Optional string, \code{"LGCM"} or \code{"LCSM"}.
#' @keywords internal
#' @noRd
validate_mgm_names <- function(y_var, y_model = NULL) {
  y_var <- as.character(y_var)
  if (anyDuplicated(y_var)) stop("The process prefixes in 'y_var' must be distinct.", call. = FALSE)
  if (identical(y_model, "LCSM") && anyDuplicated(tolower(y_var))) {
    stop("The process prefixes in 'y_var' must differ in more than letter case for a latent change score model.",
         call. = FALSE)
  }
  if (length(y_var) < 2L) return(invisible(TRUE))
  idx <- utils::combn(seq_along(y_var), 2)
  keys <- apply(idx, 2, function(i) paste0(y_var[i], collapse = ""))
  pairs <- apply(idx, 2, function(i) paste0("(\"", y_var[i[1]], "\", \"", y_var[i[2]], "\")"))
  for (j in seq_along(keys)) {
    if (keys[j] %in% y_var) {
      stop("Ambiguous multivariate names: process \"", keys[j], "\" and pair ", pairs[j], " both use the identifier \"",
           keys[j], "\". Rename one process prefix.", call. = FALSE)
    }
  }
  dup <- which(duplicated(keys))
  if (length(dup)) {
    first <- match(keys[dup[1]], keys)
    stop("Ambiguous multivariate names: pairs ", pairs[first], " and ", pairs[dup[1]], " both use the identifier \"",
         keys[dup[1]], "\". Rename one process prefix.", call. = FALSE)
  }
  invisible(TRUE)
}
