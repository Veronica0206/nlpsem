test_that(".getGF_info returns correct structure for LIN", {
  info <- nlpsem:::.getGF_info("LIN", intrinsic = FALSE)
  expect_equal(info$base_labels, c("0", "1"))
  expect_equal(info$psi_labels, c("0", "1"))
  expect_null(info$gamma_name)
})

test_that(".getGF_info returns correct structure for EXP intrinsic", {
  info <- nlpsem:::.getGF_info("EXP", intrinsic = TRUE)
  expect_equal(info$base_labels, c("0", "1"))
  expect_equal(info$psi_labels, c("0", "1", "g"))
  expect_equal(info$gamma_name, "slp_ratio")
})

test_that(".getGF_info returns correct structure for BLS intrinsic", {
  info <- nlpsem:::.getGF_info("BLS", intrinsic = TRUE)
  expect_equal(info$base_labels, c("0", "1", "2"))
  expect_equal(info$psi_labels, c("0", "1", "2", "g"))
  expect_equal(info$gamma_name, "knot")
})

test_that(".getGF_info errors on unsupported curveFun", {
  expect_error(nlpsem:::.getGF_info("cubic"), "Unsupported curveFun")
})

test_that(".lower_tri_labels generates correct pairs", {
  result <- nlpsem:::.lower_tri_labels(c("0", "1"))
  expect_equal(result, c("00", "01", "11"))

  result3 <- nlpsem:::.lower_tri_labels(c("0", "1", "2"))
  expect_equal(result3, c("00", "01", "02", "11", "12", "22"))
})

test_that(".upper_tri_labels generates correct pairs", {
  result <- nlpsem:::.upper_tri_labels(c("0", "1"))
  expect_equal(result, c("00", "01", "11"))
})

test_that(".auto_names_LGCM returns correct names for linear no TIC", {
  nms <- nlpsem:::.auto_names_LGCM("LIN", intrinsic = FALSE, growth_TIC = NULL)
  expect_true("mueta0" %in% nms)
  expect_true("mueta1" %in% nms)
  expect_true("psi00" %in% nms)
  expect_true("residuals" %in% nms)
  # mueta0, mueta1, psi00, psi01, psi11, residuals = 6
  expect_equal(length(nms), 6)
})

test_that(".auto_names_LGCM returns correct names for BLS with TIC", {
  nms <- nlpsem:::.auto_names_LGCM("BLS", intrinsic = TRUE, growth_TIC = c("ex1"))
  expect_true("alpha0" %in% nms)
  expect_true("alphag" %in% nms)
  expect_true("beta10" %in% nms)
  expect_true("mux1" %in% nms)
  expect_true("residuals" %in% nms)
})

test_that(".auto_names_subModel dispatches correctly", {
  nms_lgcm <- nlpsem:::.auto_names_subModel("LGCM", "LIN", FALSE, NULL, 1:9, NULL, NULL, "M", NULL, NULL, NULL)
  expect_true("mueta0" %in% nms_lgcm)

  expect_error(nlpsem:::.auto_names_subModel("UNKNOWN", "LIN", FALSE, NULL, 1:9, NULL, NULL, "M", NULL, NULL, NULL),
               "Unsupported sub_Model")
})
