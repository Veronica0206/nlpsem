#' @title Container-specific components for class-specific models
#' @description The sub-model builders produce one \code{mxModel} per class; the container decides what each class
#' model carries in addition. A multiple-group model gives class k its own rows of the data; a mixture model gives
#' every class the row-wise fit function and attaches the data once at the top level.
#' @param dat The data set.
#' @param nClass Number of classes.
#' @param container \code{"multigroup"} or \code{"mixture"}.
#' @param grp_var Grouping variable (codes 1..nClass), required for \code{"multigroup"}.
#' @return A list of length \code{nClass}; element k is the list of extra components for class k.
#' @keywords internal
#' @noRd
#'
#' @importFrom OpenMx mxData mxFitFunctionML
#'
.class_extras <- function(dat, nClass, container = c("multigroup", "mixture"), grp_var = NULL){
  container <- match.arg(container)
  if (container == "multigroup"){
    if (is.null(grp_var)) stop("grp_var is required for a multiple-group model.")
    lapply(seq_len(nClass), function(k) list(mxData(observed = dat[dat[, grp_var] == k, ], type = "raw")))
  } else{
    lapply(seq_len(nClass), function(k) list(mxFitFunctionML(vector = TRUE)))
  }
}

## Helper: the name a single-group model carries (kept identical to the names the 0.4 builders used).
.curve_label <- function(curveFun){
  # a lookup vector rather than switch(): an alternative named EXP partially matches switch's EXPR argument
  labels <- c("linear" = "Linear", "LIN" = "Linear", "quadratic" = "Quadratic", "QUAD" = "Quadratic",
              "negative exponential" = "Negative Exponential", "EXP" = "Negative Exponential",
              "Jenss-Bayley" = "Jenss Bayley", "JB" = "Jenss Bayley",
              "bilinear spline" = "Bilinear Spline", "BLS" = "Bilinear Spline",
              "nonparametric" = "Nonparametric", "NonP" = "Nonparametric")
  if (!curveFun %in% names(labels)) stop("unknown curveFun: ", curveFun)
  unname(labels[[curveFun]])
}
.single_model_name <- function(sub_Model, curveFun, intrinsic = FALSE, growth_TIC = NULL){
  shape <- if (curveFun %in% c("negative exponential", "EXP", "Jenss-Bayley", "JB")){
    paste0(", ", if (isTRUE(intrinsic)) "Random" else "Fixed", " Ratio")
  } else if (curveFun %in% c("bilinear spline", "BLS")){
    paste0(", ", if (isTRUE(intrinsic)) "Random" else "Fixed", " knot")
  } else ""
  paste0(sub_Model, if (!is.null(growth_TIC)) "_TICs" else "", ", ", .curve_label(curveFun), " Functional Form", shape)
}

## Helper: names of the single-group TVC, multivariate and mediation models (systematic form of the
## names the 0.4 builders used; the few copy-paste variants of 0.4 are not reproduced).
.single_model_name_tvc <- function(y_model, curveFun, intrinsic, growth_TIC, decompose){
  shape <- if (curveFun %in% c("negative exponential", "EXP", "Jenss-Bayley", "JB")){
    paste0(" (", if (isTRUE(intrinsic)) "random" else "fixed", " ratio)")
  } else if (curveFun %in% c("bilinear spline", "BLS")){
    paste0(" (", if (isTRUE(intrinsic)) "random" else "fixed", " knot)")
  } else ""
  paste0(y_model, " with ", if (!is.null(decompose) && decompose != 0) "decomposed TVC" else "TVC",
         if (!is.null(growth_TIC)) " and TICs" else "", ", ", .curve_label(curveFun), " Functional Form", shape)
}
.single_model_name_mgm <- function(y_model, curveFun, intrinsic){
  shape <- if (curveFun %in% c("negative exponential", "EXP", "Jenss-Bayley", "JB")){
    paste0(" (", if (isTRUE(intrinsic)) "random" else "fixed", " ratio)")
  } else if (curveFun %in% c("bilinear spline", "BLS")){
    paste0(" (", if (isTRUE(intrinsic)) "random" else "fixed", " knot)")
  } else ""
  paste0("Multivariate Latent ", if (y_model == "LGCM") "Growth Curve" else "Change Score", " Model, ",
         .curve_label(curveFun), shape)
}
.single_model_name_med <- function(x_type, curveFun){
  paste0("Longitudinal Mediation Model with ", if (x_type == "baseline") "Baseline" else "Longitudinal",
         " Predictor, ", .curve_label(curveFun))
}
