#' Validate that parameter names are provided when paramOut is TRUE
#' @param paramOut Logical flag for parameter output.
#' @param names Character vector of parameter names.
#' @keywords internal
#' @noRd
validate_paramOut <- function(paramOut, names) {
  if (!is.logical(paramOut) || length(paramOut) != 1) {
    stop("paramOut must be a single logical value (TRUE or FALSE).")
  }
  if (!is.null(names) && !is.character(names)) {
    stop("names must be a character vector or NULL.")
  }
  invisible(NULL)
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
