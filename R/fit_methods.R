#' Standard methods for fitted nlpsem models
#'
#' @description Standard R methods for the \code{myMxOutput} objects returned by
#' the nlpsem fitting functions. The existing S4 slots are unchanged, so saved
#' objects from earlier versions can still be summarized and inspected.
#'
#' @param x,object A fitted \code{myMxOutput} object.
#' @param digits Number of significant digits used by \code{print()}.
#' @param y Unused. Supply plotting options by name.
#' @param ... Arguments forwarded to the underlying OpenMx \code{summary()}
#' method. For \code{plot()} on a legacy object without saved specifications,
#' named structural arguments accepted by \code{getFigure()}, such as
#' \code{sub_Model}, \code{y_var}, \code{curveFun}, \code{t_var}, and \code{records}.
#' Otherwise unused by the inspection methods. Structural plotting overrides
#' are rejected for objects with a saved specification.
#' @param xstarts A numeric shift of the plotted time axis. Default is zero;
#' the package does not infer time transformations made before fitting.
#' @param xlab Label for the time axis. Default is \code{"Time"}.
#' @param outcome Display names for the longitudinal processes. By default,
#' uses the saved outcome, mediator, and longitudinal predictor prefixes.
#' @param draw Whether to display the generated figures. Set \code{FALSE} to
#' obtain the figure object without drawing it.
#'
#' @details \code{coef()} returns the named, unrounded free parameters from
#' the underlying OpenMx model, including parameters that are not shown in the
#' derived \code{Estimates} table. \code{logLik()} delegates to OpenMx and retains
#' its parameter-count and observation-count attributes for \code{AIC()} and
#' \code{BIC()}. \code{summary()} delegates to OpenMx without refitting the model.
#'
#' \code{plot()} uses \code{getSpecification()} to recover the fitted model's
#' structural arguments and forwards them to \code{getFigure()}. Only display
#' arguments (\code{xstarts}, \code{xlab}, \code{outcome}, and \code{draw}) may be
#' overridden. Legacy fits require explicit plotting specifications; the method
#' never attempts to reconstruct a missing specification from the model name.
#' Original group labels remain available in \code{getSpecification()}; the
#' concise print method reports the mapping to the fitted group codes.
#'
#' @return \code{print()} and \code{show()} return the fitted object invisibly.
#' \code{summary()} returns the underlying OpenMx summary object.
#' \code{coef()} returns a named numeric vector of free parameter estimates.
#' \code{logLik()} returns an object of class \code{logLik}, with \code{df}
#' (number of free parameters) and \code{nobs} (number of observations when
#' available) attributes. \code{plot()} returns a \code{figOutput} object invisibly,
#' drawing the figures first when \code{draw = TRUE}.
#' @seealso \code{\link{getSpecification}}, \code{\link{getFigure}},
#' \code{\link{ModelSummary}}, \code{\link{printTable}}
#' @name fit-methods
#' @examples
#' load(system.file("extdata", "getLGCM_examples.RData", package = "nlpsem"))
#' print(Math_LGCM_BLS_r)
#' head(coef(Math_LGCM_BLS_r))
#' logLik(Math_LGCM_BLS_r)
#' # The fit records its specification, so plotting takes only presentation choices:
#' figure <- plot(Math_LGCM_BLS_r, xlab = "Month", outcome = "Mathematics", draw = FALSE)
NULL

#' @rdname fit-methods
#' @export
print.myMxOutput <- function(x, digits = 4, ...){
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      !is.finite(digits) || digits < 1){
    stop("'digits' must be a positive finite number.", call. = FALSE)
  }
  model <- x@mxOutput
  information <- getSpecification(x)
  cat("nlpsem fit: ", model@name, "\n", sep = "")
  if (!is.null(information)){
    cat("Fitted with: ", information$function_name, "()\n", sep = "")
    mapping <- information$data$group_mapping
    if (is.data.frame(mapping) && nrow(mapping) > 0L && all(c("code", "label") %in% names(mapping))){
      cat("Groups: ", paste(paste0(mapping$code, " = ", mapping$label), collapse = "; "), "\n", sep = "")
    }
  }
  likelihood <- tryCatch(stats::logLik(model), error = function(e) NULL)
  n <- if (is.null(likelihood)) NULL else attr(likelihood, "nobs")
  cat("Observations: ", if (length(n) == 1L && is.finite(n)) format(n) else "unavailable", "\n", sep = "")
  cat("Free parameters: ", length(coef.myMxOutput(x)), "\n", sep = "")
  fit <- if (is.null(likelihood)) NA_real_ else as.numeric(likelihood)
  if (length(fit) == 1L && is.finite(fit)){
    cat("-2 log-likelihood: ", format(-2 * fit, digits = digits), "\n", sep = "")
  }
  status <- model@output$status$code
  cat("Optimizer status: ", if (length(status) == 1L) status else "unavailable", "\n", sep = "")
  invisible(x)
}

#' @rdname fit-methods
#' @export
summary.myMxOutput <- function(object, ...){
  summary(object@mxOutput, ...)
}

#' @rdname fit-methods
#' @export
#' @importFrom stats coef logLik
coef.myMxOutput <- function(object, ...){
  OpenMx::omxGetParameters(object@mxOutput, free = TRUE)
}

#' @rdname fit-methods
#' @export
logLik.myMxOutput <- function(object, ...){
  stats::logLik(object@mxOutput)
}

#' @rdname fit-methods
#' @export
#' @importFrom graphics plot
plot.myMxOutput <- function(x, y, ..., xstarts = 0, xlab = "Time", outcome = NULL, draw = TRUE){
  if (!missing(y)) stop("Supply plot display options by name; 'y' is not used for an nlpsem fit.", call. = FALSE)
  if (!is.logical(draw) || length(draw) != 1L || is.na(draw)){
    stop("'draw' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(xstarts) || length(xstarts) != 1L || !is.finite(xstarts)){
    stop("'xstarts' must be one finite numeric value.", call. = FALSE)
  }
  if (!is.character(xlab) || length(xlab) != 1L || is.na(xlab)){
    stop("'xlab' must be one character string.", call. = FALSE)
  }
  supplied <- list(...)
  if (length(supplied) && (is.null(names(supplied)) || any(!nzchar(names(supplied))) || anyDuplicated(names(supplied)))){
    stop("Additional plotting arguments must have unique names.", call. = FALSE)
  }
  figure_fields <- setdiff(names(formals(getFigure)), c("model", "xstarts", "xlab", "outcome"))
  information <- getSpecification(x)
  if (is.null(information)){
    unknown <- setdiff(names(supplied), figure_fields)
    if (length(unknown)) stop("Unknown plotting argument(s): ", paste(unknown, collapse = ", "), ".", call. = FALSE)
    specification <- supplied
  } else {
    if (length(supplied)){
      stop("Structural arguments are taken from the saved fit specification; remove: ",
           paste(names(supplied), collapse = ", "), ". Only xstarts, xlab, outcome, and draw may be overridden.", call. = FALSE)
    }
    specification <- information$specification
  }
  required <- c("sub_Model", "y_var", "curveFun", "t_var", "records")
  if (!is.null(specification$sub_Model) && specification$sub_Model %in% c("TVC", "MGM")) required <- c(required, "y_model")
  if (identical(specification$sub_Model, "MED")) required <- c(required, "m_var", "x_type", "x_var")
  absent <- required[vapply(required, function(nm) is.null(specification[[nm]]), logical(1))]
  if (length(absent)){
    if (is.null(information)){
      stop("This legacy fit has no saved specification. Supply the getFigure() arguments explicitly: ",
           paste(absent, collapse = ", "), ".", call. = FALSE)
    }
    stop("The saved specification is incomplete for plotting: ", paste(absent, collapse = ", "), ".", call. = FALSE)
  }
  if (is.null(outcome)){
    outcome <- specification$y_var
    if (identical(specification$sub_Model, "MED")){
      outcome <- c(outcome, specification$m_var)
      if (identical(specification$x_type, "longitudinal")) outcome <- c(outcome, specification$x_var)
    }
  }
  if (!is.character(outcome) || !length(outcome) || anyNA(outcome)){
    stop("'outcome' must contain non-missing character display names.", call. = FALSE)
  }
  args <- specification[intersect(figure_fields, names(specification))]
  figure <- do.call(getFigure, c(list(model = x@mxOutput), args,
                                 list(xstarts = xstarts, xlab = xlab, outcome = outcome)))
  if (draw) methods::show(figure)
  invisible(figure)
}
