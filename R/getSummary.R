#' @title Summarize Model Fit Statistics for Fitted Models
#'
#' @description This function summarizes the model fit statistics for a list of fitted models. The summary includes the
#' number of parameters, estimated likelihood (-2ll), AIC, BIC, and other relevant statistics.
#'
#' @param model_list A list of fitted mxModel objects.
#' @param HetModels A logical flag indicating whether a mixture model or a multiple group model is included in the list.
#' If set to \code{TRUE}, the function can also be used for the enumeration process, allowing the determination of the
#' optimal number of latent classes based on model fit statistics such as BIC. The default value is \code{FALSE}.
#'
#' @return A data frame summarizing model fit statistics (number of parameters, estimated likelihood, AIC, BIC, etc.)
#' for each model.
#'
#' @export
#'
#' @examples
#' OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
#' # Load ECLS-K (2011) data
#' data("RMS_dat")
#' RMS_dat0 <- RMS_dat
#' # Re-baseline the data so that the estimated initial status is for the starting point of the study
#' baseT <- RMS_dat0$T1
#' RMS_dat0$T1 <- RMS_dat0$T1 - baseT
#' RMS_dat0$T2 <- RMS_dat0$T2 - baseT
#' RMS_dat0$T3 <- RMS_dat0$T3 - baseT
#' RMS_dat0$T4 <- RMS_dat0$T4 - baseT
#' RMS_dat0$T5 <- RMS_dat0$T5 - baseT
#' RMS_dat0$T6 <- RMS_dat0$T6 - baseT
#' RMS_dat0$T7 <- RMS_dat0$T7 - baseT
#' RMS_dat0$T8 <- RMS_dat0$T8 - baseT
#' RMS_dat0$T9 <- RMS_dat0$T9 - baseT
#' # Fit bilinear spline growth model with fix knot (non-intrinsically nonlinear model)
#' \donttest{
#' BLS_LGCM1 <- getLGCM(
#'   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = FALSE,
#'   records = 1:9, res_scale = 0.1
#'   )
#' getSummary(model_list = list(BLS_LGCM1), HetModels = FALSE)
#' set.seed(20191029)
#' BLS_LGCM2 <-  getMIX(
#'   dat = RMS_dat0, prop_starts = c(0.45, 0.55), sub_Model = "LGCM", cluster_TIC = NULL,
#'   y_var = "M", t_var = "T", records = 1:9, curveFun = "BLS", intrinsic = FALSE,
#'   res_scale = list(0.3, 0.3), growth_TIC = NULL, tries = 10
#'   )
#' set.seed(20191029)
#' BLS_LGCM3 <-  getMIX(
#'   dat = RMS_dat0, prop_starts = c(0.33, 0.34, 0.33), sub_Model = "LGCM", cluster_TIC = NULL,
#'   y_var = "M", t_var = "T", records = 1:9, curveFun = "BLS", intrinsic = FALSE,
#'   res_scale = list(0.3, 0.3, 0.3), growth_TIC = NULL, tries = 10
#'   )
#'
#' getSummary(model_list = list(BLS_LGCM1, BLS_LGCM2, BLS_LGCM3), HetModels = TRUE)
#' }
#'
#' @importFrom dplyr left_join
#' @importFrom tidyr separate spread pivot_wider
#' @importFrom stringr str_replace_all
#' @importFrom stats AIC BIC
#'
getSummary <- function(model_list, HetModels = FALSE){
  if (HetModels){
    output0 <- list()
    for (model in 1:length(model_list)){
      params <- omxGetParameters(model_list[[model]])
      res_indices <- grep("_res|_RES", names(params))
      res_names <- names(params)[res_indices]
      res_estimates <- round(params[res_indices], 4)
      nClass <- length(model_list[[model]]@submodels)
      if (nClass == 0){
        output0[[model]] <- data.frame(Model = paste0("Model", model),
                                       No_Params = length(params),
                                       `-2ll` = round(model_list[[model]]@output$Minus2LogLikelihood, 2),
                                       AIC = round(AIC(model_list[[model]]), 2),
                                       BIC = round(BIC(model_list[[model]]), 2),
                                       Class = "c1", Prop = "100%")
        output0[[model]][1, res_names] <- res_estimates
      }
      else if (nClass > 0){
        if (is.null(model_list[[model]]$data$observed)){
          model_names <- names(model_list[[model]]@submodels)
          sub_dat <- lapply(model_names, function(name) {
            as.data.frame(model_list[[model]][[name]]$data$observed)
          })
          sub_cnt <- lapply(sub_dat, function(dat) nrow(dat))
          cnt <- unlist(sub_cnt)
          prop <- paste0(cnt/sum(cnt) * 100, "%")
        }
        else{
          cluster_TIC <- if (is.null(model_list[[model]]$weightsV)){
            NULL
          }else{
            substring(model_list[[model]]$weightsV$labels[-1], 6)
          }
          memebership <- getPosterior(model = model_list[[model]], nClass = nClass, cluster_TIC = cluster_TIC)$membership
          prop <- paste0(table(memebership)/length(memebership) * 100, "%")
        }
        output0[[model]] <- data.frame(Model = paste0("Model", model),
                                       No_Params = length(params),
                                       `-2ll` = round(model_list[[model]]@output$Minus2LogLikelihood, 2),
                                       AIC = round(AIC(model_list[[model]]), 2),
                                       BIC = round(BIC(model_list[[model]]), 2),
                                       Class = paste0("c", 1:nClass), Prop = prop)
        tmp_l <- data.frame(res_names, res_estimates)
        tmp_l <- separate(tmp_l, col = "res_names", into = c("Class", "Residual"), sep = "(?<=^.{2})")
        tmp_w <- spread(tmp_l, key = "Residual", value = "res_estimates")
        output0[[model]] <- output0[[model]] %>% left_join(tmp_w, by = "Class")
      }
    }
    output1 <- do.call(rbind, output0)
    res_cols <- grep("_residuals$", names(output1), value = TRUE)
    output_list <- lapply(res_cols, function(res_cols) {
      output1 %>%
        select(Model, No_Params, X.2ll, AIC, BIC, Class, Prop, matches(res_cols)) %>%
        pivot_wider(names_from = Class,
                    values_from = c(matches(res_cols), Prop),
                    names_sep = "_")
    })
    output <- Reduce(function(x, y) {
      merge(x, y, by = c("Model", "No_Params", "X.2ll", "AIC", "BIC"))
    }, output_list)
    names(output) <- str_replace_all(names(output), "_residuals_", "_res_")
    names(output) <- str_replace_all(names(output), "Prop_c", "%Class")
  }
  else if (!HetModels){
    # create an empty data frame with columns for the output
    output0 <- list()
    # loop over the models and extract the relevant information
    for (model in 1:length(model_list)){
      params <- omxGetParameters(model_list[[model]])
      res_indices <- grep("_res|_RES", names(params))
      res_names <- names(params)[res_indices]
      res_estimates <- round(params[res_indices], 4)
      output0[[model]] <- data.frame(Model = paste0("Model", model),
                                     No_Params = length(params),
                                     `-2ll` = round(model_list[[model]]@output$Minus2LogLikelihood, 2),
                                     AIC = round(AIC(model_list[[model]]), 2),
                                     BIC = round(BIC(model_list[[model]]), 2))
      output0[[model]][1, res_names] <- res_estimates
    }
    output <- do.call(rbind, output0)
  }
  names(output)[3] <- "-2ll"
  return(output)
}
