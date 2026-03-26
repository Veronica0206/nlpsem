# Auto-Generate Parameter Names for Model Output
#
# Internal helper functions to automatically generate meaningful parameter names
# when the user does not supply a `names` argument with `paramOut = TRUE`.

# ---- Core helpers ----

#' Get growth factor structural info for a given functional form
#' @return A list with base_labels (non-gamma GFs), psi_labels (for psi matrix),
#'   and gamma_name (semantic name for the nonlinear parameter, or NULL)
#' @keywords internal
#' @noRd
.getGF_info <- function(curveFun, intrinsic = TRUE) {
  if (curveFun %in% c("linear", "LIN")) {
    list(base_labels = c("0", "1"), psi_labels = c("0", "1"), gamma_name = NULL)
  } else if (curveFun %in% c("quadratic", "QUAD")) {
    list(base_labels = c("0", "1", "2"), psi_labels = c("0", "1", "2"), gamma_name = NULL)
  } else if (curveFun %in% c("negative exponential", "EXP")) {
    if (intrinsic) {
      list(base_labels = c("0", "1"), psi_labels = c("0", "1", "g"), gamma_name = "slp_ratio")
    } else {
      list(base_labels = c("0", "1"), psi_labels = c("0", "1"), gamma_name = "slp_ratio")
    }
  } else if (curveFun %in% c("Jenss-Bayley", "JB")) {
    if (intrinsic) {
      list(base_labels = c("0", "1", "2"), psi_labels = c("0", "1", "2", "g"), gamma_name = "acc_ratio")
    } else {
      list(base_labels = c("0", "1", "2"), psi_labels = c("0", "1", "2"), gamma_name = "acc_ratio")
    }
  } else if (curveFun %in% c("bilinear spline", "BLS")) {
    if (intrinsic) {
      list(base_labels = c("0", "1", "2"), psi_labels = c("0", "1", "2", "g"), gamma_name = "knot")
    } else {
      list(base_labels = c("0", "1", "2"), psi_labels = c("0", "1", "2"), gamma_name = "knot")
    }
  } else if (curveFun %in% c("nonparametric", "NonP")) {
    list(base_labels = c("0", "1"), psi_labels = c("0", "1"), gamma_name = NULL)
  } else {
    stop("Unsupported curveFun: ", curveFun)
  }
}

#' Generate lower-triangular psi/phi label pairs (column-major)
#' @param labels Character vector of index labels
#' @return Character vector of paired labels
#' @keywords internal
#' @noRd
.lower_tri_labels <- function(labels) {
  n <- length(labels)
  result <- c()
  for (j in seq_len(n)) {
    for (i in j:n) {
      result <- c(result, paste0(labels[j], labels[i]))
    }
  }
  result
}

#' Generate upper-triangular label pairs (row-major)
#' @param labels Character vector of index labels
#' @return Character vector of paired labels
#' @keywords internal
#' @noRd
.upper_tri_labels <- function(labels) {
  n <- length(labels)
  result <- c()
  for (i in seq_len(n)) {
    for (j in i:n) {
      result <- c(result, paste0(labels[i], labels[j]))
    }
  }
  result
}

# ---- Gamma parameter label helpers ----

#' Get the gamma mean label for no-TIC context (standalone LGCM/LCSM)
#' Convention: EXP-F = "slp_ratio", EXP-T = "mu_slp_ratio",
#'   JB-F/T = "acc_ratio" (never mu_), BLS-F = "knot", BLS-T = "mu_knot"
#' @keywords internal
#' @noRd
.gamma_noTIC_label <- function(curveFun, intrinsic) {
  if (curveFun %in% c("negative exponential", "EXP")) {
    if (intrinsic) "mu_slp_ratio" else "slp_ratio"
  } else if (curveFun %in% c("Jenss-Bayley", "JB")) {
    "acc_ratio"
  } else if (curveFun %in% c("bilinear spline", "BLS")) {
    if (intrinsic) "mu_knot" else "knot"
  } else {
    NULL
  }
}

#' Get the gamma label in alpha position (with-TIC or decompose!=0)
#' Convention: intrinsic=F uses specific name, intrinsic=T uses "alphag"
#' @keywords internal
#' @noRd
.gamma_alpha_label <- function(curveFun, intrinsic) {
  gn <- .gamma_noTIC_label(curveFun, FALSE)  # base name without mu_
  if (is.null(gn)) return(NULL)
  if (intrinsic) "alphag" else gn
}

#' Get the gamma label in mueta position (with-TIC, intrinsic=T only)
#' Convention: always "mu_" + gamma_name
#' @keywords internal
#' @noRd
.gamma_mueta_label <- function(curveFun, intrinsic) {
  if (!intrinsic) return(NULL)
  gf <- .getGF_info(curveFun, intrinsic)
  if (is.null(gf$gamma_name)) return(NULL)
  paste0("mu_", gf$gamma_name)
}

#' Get the gamma mean label for TVC decompose=0 no-TIC context
#' Convention: EXP always uses "mu_slp_ratio", JB uses "acc_ratio",
#'   BLS-F = "knot", BLS-T = "mu_knot"
#' @keywords internal
#' @noRd
.gamma_noTIC_label_tvc <- function(curveFun, intrinsic) {
  if (curveFun %in% c("negative exponential", "EXP")) {
    "mu_slp_ratio"
  } else {
    .gamma_noTIC_label(curveFun, intrinsic)
  }
}

# ---- LGCM names ----

#' Auto-generate parameter names for getLGCM
#' @keywords internal
#' @noRd
.auto_names_LGCM <- function(curveFun, intrinsic, growth_TIC) {
  gf <- .getGF_info(curveFun, intrinsic)
  psi_ltri <- .lower_tri_labels(gf$psi_labels)
  nTIC <- if (!is.null(growth_TIC)) length(growth_TIC) else 0

  if (nTIC > 0) {
    # Alpha labels
    alpha <- paste0("alpha", gf$base_labels)
    gamma_a <- .gamma_alpha_label(curveFun, intrinsic)
    if (!is.null(gamma_a)) alpha <- c(alpha, gamma_a)

    # Beta labels: beta{TIC_idx}{GF_idx} for psi (random) GFs
    beta <- unlist(lapply(seq_len(nTIC), function(t) paste0("beta", t, gf$psi_labels)))

    # TIC parameters
    phi_ltri <- .lower_tri_labels(as.character(seq_len(nTIC)))
    mux <- paste0("mux", seq_len(nTIC))
    phi <- paste0("phi", phi_ltri)

    # Mueta labels (random GFs only)
    mueta <- paste0("mueta", gf$base_labels)
    gamma_m <- .gamma_mueta_label(curveFun, intrinsic)
    if (!is.null(gamma_m)) mueta <- c(mueta, gamma_m)

    c(alpha, paste0("psi", psi_ltri), "residuals", beta, mux, phi, mueta)
  } else {
    # Mean labels
    mean_lab <- paste0("mueta", gf$base_labels)
    gamma_lab <- .gamma_noTIC_label(curveFun, intrinsic)
    if (!is.null(gamma_lab)) mean_lab <- c(mean_lab, gamma_lab)

    c(mean_lab, paste0("psi", psi_ltri), "residuals")
  }
}

# ---- LCSM names ----

#' Auto-generate parameter names for getLCSM
#' @keywords internal
#' @noRd
.auto_names_LCSM <- function(curveFun, intrinsic, growth_TIC, nT) {
  gf <- .getGF_info(curveFun, intrinsic)
  psi_ltri <- .lower_tri_labels(gf$psi_labels)
  nTIC <- if (!is.null(growth_TIC)) length(growth_TIC) else 0
  nIntervals <- nT - 1
  is_nonp <- curveFun %in% c("nonparametric", "NonP")

  if (nTIC > 0) {
    alpha <- paste0("alpha", gf$base_labels)
    gamma_a <- .gamma_alpha_label(curveFun, intrinsic)
    if (!is.null(gamma_a)) alpha <- c(alpha, gamma_a)

    rel_rate <- if (is_nonp) paste0("rel_rate", 2:(nT - 1)) else NULL

    beta <- unlist(lapply(seq_len(nTIC), function(t) paste0("beta", t, gf$psi_labels)))
    phi_ltri <- .lower_tri_labels(as.character(seq_len(nTIC)))
    mux <- paste0("mux", seq_len(nTIC))
    phi <- paste0("phi", phi_ltri)

    mueta <- paste0("mueta", gf$base_labels)
    gamma_m <- .gamma_mueta_label(curveFun, intrinsic)
    if (!is.null(gamma_m)) mueta <- c(mueta, gamma_m)

    core <- c(alpha, paste0("psi", psi_ltri), rel_rate, "residuals",
              beta, mux, phi, mueta)
  } else {
    mean_lab <- paste0("mueta", gf$base_labels)
    gamma_lab <- .gamma_noTIC_label(curveFun, intrinsic)
    if (!is.null(gamma_lab)) mean_lab <- c(mean_lab, gamma_lab)

    rel_rate <- if (is_nonp) paste0("rel_rate", 2:(nT - 1)) else NULL

    core <- c(mean_lab, paste0("psi", psi_ltri), rel_rate, "residuals")
  }

  # Append derived change-score quantities (no prefix for standalone LCSM)
  derived <- c(paste0("slp_val_est", seq_len(nIntervals)),
               paste0("slp_var_est", seq_len(nIntervals)),
               paste0("chg_inv_val_est", seq_len(nIntervals)),
               paste0("chg_inv_var_est", seq_len(nIntervals)),
               paste0("chg_bl_val_est", seq_len(nIntervals)),
               paste0("chg_bl_var_est", seq_len(nIntervals)))

  c(core, derived)
}

# ---- TVC names ----

#' Auto-generate parameter names for getTVCmodel
#' @keywords internal
#' @noRd
.auto_names_TVC <- function(curveFun, intrinsic, records, y_model, decompose, growth_TIC) {
  gf <- .getGF_info(curveFun, intrinsic)
  psi_ltri <- .lower_tri_labels(gf$psi_labels)
  nTIC <- if (!is.null(growth_TIC)) length(growth_TIC) else 0
  nT <- length(records)
  nIntervals <- nT - 1
  is_nonp <- curveFun %in% c("nonparametric", "NonP")

  if (decompose == 0) {
    .auto_names_TVC_d0(gf, curveFun, intrinsic, psi_ltri, records, nT,
                        nIntervals, y_model, nTIC, is_nonp)
  } else {
    .auto_names_TVC_dx(gf, curveFun, intrinsic, psi_ltri, records, nT,
                        nIntervals, y_model, nTIC, is_nonp)
  }
}

#' TVC decompose=0 names
#' @keywords internal
#' @noRd
.auto_names_TVC_d0 <- function(gf, curveFun, intrinsic, psi_ltri, records, nT,
                                 nIntervals, y_model, nTIC, is_nonp) {
  # Y-model parameters (with Y_ prefix)
  if (nTIC > 0) {
    y_alpha <- paste0("Y_alpha", gf$base_labels)
    gamma_a <- .gamma_alpha_label(curveFun, intrinsic)
    if (!is.null(gamma_a)) y_alpha <- c(y_alpha, paste0("Y_", gamma_a))
    y_mean <- y_alpha
  } else {
    y_mean <- paste0("Y_mueta", gf$base_labels)
    gamma_lab <- .gamma_noTIC_label_tvc(curveFun, intrinsic)
    if (!is.null(gamma_lab)) y_mean <- c(y_mean, paste0("Y_", gamma_lab))
  }

  y_psi <- paste0("Y_psi", psi_ltri)
  rel_rate <- if (is_nonp && y_model == "LCSM") paste0("Y_rel_rate", 2:(nT - 1)) else NULL

  # TVC parameters
  tvc_m <- paste0("TVC_m", records)
  tvc_v <- paste0("TVC_v", records)

  # TIC parameters (if any)
  tic_params <- NULL
  if (nTIC > 0) {
    if (nTIC == 1) {
      tic_params <- c(paste0("betaTIC", gf$psi_labels),
                       "muTIC", "phiTIC")
    } else {
      beta_tic <- unlist(lapply(seq_len(nTIC), function(t) paste0("betaTIC", t, gf$psi_labels)))
      phi_ltri <- .lower_tri_labels(as.character(seq_len(nTIC)))
      tic_params <- c(beta_tic,
                       paste0("muTIC", seq_len(nTIC)),
                       paste0("phiTIC", phi_ltri))
    }
    # Mueta section
    y_mueta <- paste0("Y_mueta", gf$base_labels)
    gamma_m <- .gamma_mueta_label(curveFun, intrinsic)
    if (!is.null(gamma_m)) y_mueta <- c(y_mueta, paste0("Y_", gamma_m))
  } else {
    y_mueta <- NULL
  }

  if (y_model == "LGCM") {
    c(y_mean, y_psi, "Y_residuals", tvc_m, tvc_v,
      tic_params, y_mueta, "kappa")
  } else {
    # LCSM
    core <- c(y_mean, y_psi, rel_rate, "Y_residuals", tvc_m, tvc_v,
              tic_params, y_mueta, "kappa")
    derived <- c(paste0("Y_slp_val_est", seq_len(nIntervals)),
                 paste0("Y_slp_var_est", seq_len(nIntervals)),
                 paste0("Y_chg_inv_val_est", seq_len(nIntervals)),
                 paste0("Y_chg_inv_var_est", seq_len(nIntervals)),
                 paste0("Y_chg_bl_val_est", seq_len(nIntervals)),
                 paste0("Y_chg_bl_var_est", seq_len(nIntervals)))
    c(core, derived)
  }
}

#' TVC decompose!=0 names
#' @keywords internal
#' @noRd
.auto_names_TVC_dx <- function(gf, curveFun, intrinsic, psi_ltri, records, nT,
                                 nIntervals, y_model, nTIC, is_nonp) {
  # Y-model: always uses alpha convention (decompose!=0 treats like with-TIC)
  y_alpha <- paste0("Y_alpha", gf$base_labels)
  gamma_a <- .gamma_alpha_label(curveFun, intrinsic)
  if (!is.null(gamma_a)) y_alpha <- c(y_alpha, paste0("Y_", gamma_a))

  y_psi <- paste0("Y_psi", psi_ltri)
  rel_rate_y <- if (is_nonp && y_model == "LCSM") paste0("Y_rel_rate", 2:(nT - 1)) else NULL

  # X-model (NonP latent growth for TVC)
  x_psi_ltri <- .lower_tri_labels(c("0", "1"))
  x_params <- c("X_mueta0", "X_mueta1",
                 paste0("X_psi", x_psi_ltri),
                 paste0("X_rel_rate", 2:(nT - 1)),
                 paste0("X_abs_rate", records[-1]),
                 "X_residuals")

  # Beta for TVC effect
  betaTVC <- paste0("betaTVC", gf$psi_labels)

  # TIC parameters
  tic_params <- NULL
  covBL <- NULL
  if (nTIC > 0) {
    if (nTIC == 1) {
      tic_params <- c(paste0("betaTIC", gf$psi_labels), betaTVC,
                       "muTIC", "phiTIC")
      covBL <- "covBL"
    } else {
      beta_tic <- unlist(lapply(seq_len(nTIC), function(t) paste0("betaTIC", t, gf$psi_labels)))
      phi_ltri <- .lower_tri_labels(as.character(seq_len(nTIC)))
      tic_params <- c(beta_tic, betaTVC,
                       paste0("muTIC", seq_len(nTIC)),
                       paste0("phiTIC", phi_ltri))
      covBL <- paste0("covBL", seq_len(nTIC))
    }
  } else {
    tic_params <- betaTVC
    covBL <- NULL
  }

  # Mueta section
  y_mueta <- paste0("Y_mueta", gf$base_labels)
  gamma_m <- .gamma_mueta_label(curveFun, intrinsic)
  if (!is.null(gamma_m)) y_mueta <- c(y_mueta, paste0("Y_", gamma_m))

  if (y_model == "LGCM") {
    c(y_alpha, y_psi, "Y_residuals", x_params,
      tic_params, y_mueta, covBL, "kappa", "Cov_XYres")
  } else {
    # LCSM
    core <- c(y_alpha, y_psi, rel_rate_y, "Y_residuals", x_params,
              tic_params, y_mueta, covBL, "kappa", "Cov_XYres")
    derived <- c(paste0("Y_slp_val_est", seq_len(nIntervals)),
                 paste0("Y_slp_var_est", seq_len(nIntervals)),
                 paste0("Y_chg_inv_val_est", seq_len(nIntervals)),
                 paste0("Y_chg_inv_var_est", seq_len(nIntervals)),
                 paste0("Y_chg_bl_val_est", seq_len(nIntervals)),
                 paste0("Y_chg_bl_var_est", seq_len(nIntervals)))
    c(core, derived)
  }
}

# ---- MGM names ----

#' Get the gamma mean label for MGM context (per-trajectory, uses specific name without mu_)
#' @keywords internal
#' @noRd
.gamma_mgm_label <- function(curveFun) {
  if (curveFun %in% c("negative exponential", "EXP")) "slp_ratio"
  else if (curveFun %in% c("Jenss-Bayley", "JB")) "acc_ratio"
  else if (curveFun %in% c("bilinear spline", "BLS")) "knot"
  else NULL
}

#' Auto-generate parameter names for getMGM
#' @keywords internal
#' @noRd
.auto_names_MGM <- function(curveFun, intrinsic, y_var, records, y_model) {
  gf <- .getGF_info(curveFun, intrinsic)
  psi_ltri <- .lower_tri_labels(gf$psi_labels)
  nTraj <- length(y_var)
  gamma_lab <- .gamma_mgm_label(curveFun)

  # Per-trajectory parameters
  outcome_names <- c()
  for (traj in seq_len(nTraj)) {
    v <- y_var[traj]
    # Mean labels: v_mueta{base} + v_gamma
    traj_mean <- paste0(v, "_mueta", gf$base_labels)
    if (!is.null(gamma_lab)) traj_mean <- c(traj_mean, paste0(v, "_", gamma_lab))
    # Psi labels
    traj_psi <- paste0(v, "_psi", psi_ltri)
    outcome_names <- c(outcome_names, traj_mean, traj_psi, paste0(v, "_res"))
  }

  if (y_model == "LGCM") {
    # Between-trajectory cross-psi and residual
    btw_names <- .mgm_between_names(y_var, gf$psi_labels)
    c(outcome_names, btw_names)
  } else {
    # LCSM: add rel_rate, then derived (interleaved by type), then between
    is_nonp <- curveFun %in% c("nonparametric", "NonP")

    rel_names <- NULL
    if (is_nonp) {
      for (traj in seq_len(nTraj)) {
        nT_traj <- length(records[[traj]])
        rel_names <- c(rel_names, paste0(y_var[traj], "_rel_rate", 2:(nT_traj - 1)))
      }
    }

    # Derived quantities: interleaved by type across trajectories
    derived_types <- c("slp_val_est", "slp_var_est", "chg_inv_val_est",
                        "chg_inv_var_est", "chg_bl_val_est", "chg_bl_var_est")
    derived_names <- c()
    for (dtype in derived_types) {
      for (traj in seq_len(nTraj)) {
        nT_traj <- length(records[[traj]])
        nInt <- nT_traj - 1
        derived_names <- c(derived_names, paste0(y_var[traj], "_", dtype, seq_len(nInt)))
      }
    }

    btw_names <- .mgm_between_names(y_var, gf$psi_labels)

    if (is_nonp) {
      c(outcome_names, rel_names, derived_names, btw_names)
    } else {
      c(outcome_names, derived_names, btw_names)
    }
  }
}

#' Generate between-trajectory cross-psi labels for MGM
#' @keywords internal
#' @noRd
.mgm_between_names <- function(y_var, psi_labels) {
  nTraj <- length(y_var)
  nPsi <- length(psi_labels)
  btw_names <- c()
  for (i in seq_len(nTraj - 1)) {
    for (j in (i + 1):nTraj) {
      # Full cross-psi matrix (nPsi × nPsi), column-major
      cross_labels <- c()
      for (col_idx in seq_len(nPsi)) {
        for (row_idx in seq_len(nPsi)) {
          cross_labels <- c(cross_labels,
                            paste0(psi_labels[row_idx], psi_labels[col_idx]))
        }
      }
      btw_names <- c(btw_names,
                      paste0(y_var[i], y_var[j], "_psi", cross_labels),
                      paste0(y_var[i], y_var[j], "_res"))
    }
  }
  btw_names
}

# ---- MED names ----

#' Auto-generate parameter names for getMediation
#' @keywords internal
#' @noRd
.auto_names_MED <- function(curveFun, x_type, x_var, m_var, y_var) {
  if (curveFun %in% c("linear", "LIN")) {
    gf <- c("0", "1")
    if (x_type == "baseline") .med_names_LIN_baseline(gf) else .med_names_LIN_longitudinal(gf)
  } else {
    # BLS
    gf <- c("1", "r", "2")
    if (x_type == "baseline") .med_names_BLS_baseline(gf) else .med_names_BLS_longitudinal(gf)
  }
}

#' MED LIN baseline names
#' @keywords internal
#' @noRd
.med_names_LIN_baseline <- function(gf) {
  psi_ltri <- .lower_tri_labels(gf)
  upper_tri <- .upper_tri_labels(gf)
  c(
    "muX", "phi11",
    paste0("alpha", "M", gf),
    paste0("psi", paste0("M", gf[rep(seq_along(gf), each = length(gf))],
                         "M", gf[rep(seq_along(gf), length(gf))]), "_r")[
      .ltri_idx(length(gf))],
    paste0("alpha", "Y", gf),
    paste0("psi", paste0("Y", gf[rep(seq_along(gf), each = length(gf))],
                         "Y", gf[rep(seq_along(gf), length(gf))]), "_r")[
      .ltri_idx(length(gf))],
    paste0("beta", rep(c("M", "Y"), each = length(gf)), rep(gf, 2)),
    paste0("beta", .med_upper_pairs("M", "Y", gf)),
    paste0("mueta", "M", gf),
    paste0("mueta", "Y", gf),
    paste0("Mediator", upper_tri),
    paste0("total", gf),
    "residualsY", "residualsM", "residualsYM"
  )
}

#' MED LIN longitudinal names
#' @keywords internal
#' @noRd
.med_names_LIN_longitudinal <- function(gf) {
  psi_ltri <- .lower_tri_labels(gf)
  upper_tri <- .upper_tri_labels(gf)
  c(
    paste0("mueta", "X", gf),
    .med_psi_labels("X", gf, suffix = ""),
    paste0("alpha", "M", gf),
    .med_psi_labels("M", gf, suffix = "_r"),
    paste0("alpha", "Y", gf),
    .med_psi_labels("Y", gf, suffix = "_r"),
    paste0("beta", .med_upper_pairs("X", "Y", gf)),
    paste0("beta", .med_upper_pairs("X", "M", gf)),
    paste0("beta", .med_upper_pairs("M", "Y", gf)),
    paste0("mueta", "M", gf),
    paste0("mueta", "Y", gf),
    paste0("mediator", .med_3index(gf)),
    paste0("total", upper_tri),
    "residualsX", "residualsM", "residualsY", "residualsMX", "residualsYX", "residualsYM"
  )
}

#' MED BLS baseline names
#' @keywords internal
#' @noRd
.med_names_BLS_baseline <- function(gf) {
  upper_tri <- .upper_tri_labels(gf)
  c(
    "muX", "phi11",
    paste0("alpha", "M", gf), "mugM",
    .med_psi_labels("M", gf, suffix = "_r"),
    paste0("alpha", "Y", gf), "mugY",
    .med_psi_labels("Y", gf, suffix = "_r"),
    paste0("beta", rep(c("M", "Y"), each = length(gf)), rep(gf, 2)),
    paste0("beta", .med_upper_pairs("M", "Y", gf)),
    paste0("mueta", "M", gf),
    paste0("mueta", "Y", gf),
    paste0("Mediator", upper_tri),
    paste0("total", gf),
    "residualsM", "residualsY", "residualsYM"
  )
}

#' MED BLS longitudinal names
#' @keywords internal
#' @noRd
.med_names_BLS_longitudinal <- function(gf) {
  upper_tri <- .upper_tri_labels(gf)
  c(
    paste0("mueta", "X", gf), "mugX",
    .med_psi_labels("X", gf, suffix = ""),
    paste0("alpha", "M", gf), "mugM",
    .med_psi_labels("M", gf, suffix = "_r"),
    paste0("alpha", "Y", gf), "mugY",
    .med_psi_labels("Y", gf, suffix = "_r"),
    paste0("beta", .med_upper_pairs("X", "Y", gf)),
    paste0("beta", .med_upper_pairs("X", "M", gf)),
    paste0("beta", .med_upper_pairs("M", "Y", gf)),
    paste0("mueta", "M", gf),
    paste0("mueta", "Y", gf),
    paste0("mediator", .med_3index(gf)),
    paste0("total", upper_tri),
    "residualsX", "residualsM", "residualsY", "residualsMX", "residualsYX", "residualsYM"
  )
}

#' Generate variable-specific psi labels for MED (lower triangular)
#' @keywords internal
#' @noRd
.med_psi_labels <- function(var, gf, suffix = "") {
  n <- length(gf)
  result <- c()
  for (j in seq_len(n)) {
    for (i in j:n) {
      result <- c(result, paste0("psi", var, gf[j], var, gf[i], suffix))
    }
  }
  result
}

#' Generate upper-triangular beta pair labels for MED
#' @keywords internal
#' @noRd
.med_upper_pairs <- function(var1, var2, gf) {
  n <- length(gf)
  result <- c()
  for (i in seq_len(n)) {
    for (j in i:n) {
      result <- c(result, paste0(var1, gf[i], var2, gf[j]))
    }
  }
  result
}

#' Generate 3-index mediator labels for longitudinal MED
#' @keywords internal
#' @noRd
.med_3index <- function(gf) {
  n <- length(gf)
  result <- c()
  for (i in seq_len(n)) {
    for (j in i:n) {
      for (k in j:n) {
        result <- c(result, paste0(gf[i], gf[j], gf[k]))
      }
    }
  }
  result
}

#' Get lower-triangular indices for an n×n matrix
#' @keywords internal
#' @noRd
.ltri_idx <- function(n) {
  idx <- c()
  for (j in seq_len(n)) {
    for (i in j:n) {
      idx <- c(idx, (j - 1) * n + i)
    }
  }
  idx
}

# ---- Dispatchers for MIX / MGroup ----

#' Auto-generate per-class parameter names based on sub_Model
#' @keywords internal
#' @noRd
.auto_names_subModel <- function(sub_Model, curveFun, intrinsic, growth_TIC,
                                  records, y_model, decompose,
                                  y_var, m_var, x_type, x_var) {
  if (sub_Model == "LGCM") {
    .auto_names_LGCM(curveFun, intrinsic, growth_TIC)
  } else if (sub_Model == "LCSM") {
    nT <- length(records)
    .auto_names_LCSM(curveFun, intrinsic, growth_TIC, nT)
  } else if (sub_Model == "TVC") {
    .auto_names_TVC(curveFun, intrinsic, records, y_model, decompose, growth_TIC)
  } else if (sub_Model == "MGM") {
    .auto_names_MGM(curveFun, intrinsic, y_var, records, y_model)
  } else if (sub_Model == "MED") {
    .auto_names_MED(curveFun, x_type, x_var, m_var, y_var)
  } else {
    stop("Unsupported sub_Model: ", sub_Model)
  }
}
