test_that("validate_paramOut rejects non-logical paramOut", {
  expect_error(nlpsem:::validate_paramOut("TRUE", NULL), "paramOut must be a single logical")
  expect_error(nlpsem:::validate_paramOut(c(TRUE, FALSE), NULL), "paramOut must be a single logical")
  expect_error(nlpsem:::validate_paramOut(1, NULL), "paramOut must be a single logical")
})

test_that("validate_paramOut rejects non-character names", {
  expect_error(nlpsem:::validate_paramOut(TRUE, 123), "names must be a character vector")
})

test_that("validate_paramOut passes for valid inputs", {
  expect_silent(nlpsem:::validate_paramOut(TRUE, NULL))
  expect_silent(nlpsem:::validate_paramOut(FALSE, NULL))
  expect_silent(nlpsem:::validate_paramOut(TRUE, c("a", "b")))
})

test_that("validate_res_scale accepts NULL", {
  expect_silent(nlpsem:::validate_res_scale(NULL))
})

test_that("validate_res_scale rejects out-of-range values", {
  expect_error(nlpsem:::validate_res_scale(0), "between 0 and 1")
  expect_error(nlpsem:::validate_res_scale(1), "between 0 and 1")
  expect_error(nlpsem:::validate_res_scale(-0.5), "between 0 and 1")
  expect_error(nlpsem:::validate_res_scale(1.5), "between 0 and 1")
})

test_that("validate_res_scale accepts valid values", {
  expect_silent(nlpsem:::validate_res_scale(0.5))
  expect_silent(nlpsem:::validate_res_scale(c(0.1, 0.9)))
})

test_that("validate_res_scale handles list mode for getMGroup/getMIX", {
  expect_error(nlpsem:::validate_res_scale(list(0, 0.5), is_list = TRUE), "between 0 and 1")
  expect_silent(nlpsem:::validate_res_scale(list(0.1, 0.5), is_list = TRUE))
})

test_that("validate_res_cor accepts NULL", {
  expect_silent(nlpsem:::validate_res_cor(NULL))
})

test_that("validate_res_cor rejects out-of-range values", {
  expect_error(nlpsem:::validate_res_cor(-1), "between -1 and 1")
  expect_error(nlpsem:::validate_res_cor(1), "between -1 and 1")
  expect_error(nlpsem:::validate_res_cor(2), "between -1 and 1")
})

test_that("validate_res_cor accepts valid values", {
  expect_silent(nlpsem:::validate_res_cor(0.3))
  expect_silent(nlpsem:::validate_res_cor(-0.5))
  expect_silent(nlpsem:::validate_res_cor(c(0.1, -0.3)))
})

test_that("validate_intrinsic accepts NULL intrinsic", {
  expect_silent(nlpsem:::validate_intrinsic(NULL, "LIN"))
})

test_that("validate_intrinsic rejects invalid LGCM combos", {
  expect_error(nlpsem:::validate_intrinsic(TRUE, "linear", "LGCM"), "negative exponential")
  expect_error(nlpsem:::validate_intrinsic(TRUE, "LIN", "LGCM"), "negative exponential")
  expect_error(nlpsem:::validate_intrinsic(TRUE, "QUAD", "LGCM"), "negative exponential")
})

test_that("validate_intrinsic accepts valid LGCM combos", {
  expect_silent(nlpsem:::validate_intrinsic(TRUE, "EXP", "LGCM"))
  expect_silent(nlpsem:::validate_intrinsic(TRUE, "JB", "LGCM"))
  expect_silent(nlpsem:::validate_intrinsic(TRUE, "BLS", "LGCM"))
  expect_silent(nlpsem:::validate_intrinsic(FALSE, "LIN", "LGCM"))
})

test_that("validate_intrinsic rejects invalid LCSM combos", {
  expect_error(nlpsem:::validate_intrinsic(TRUE, "NonP", "LCSM"), "negative exponential")
  expect_error(nlpsem:::validate_intrinsic(TRUE, "QUAD", "LCSM"), "negative exponential")
})

test_that("validate_curveFun_ymodel rejects LCSM with unsupported forms", {
  expect_error(nlpsem:::validate_curveFun_ymodel("LIN", "LCSM"), "does not support")
  expect_error(nlpsem:::validate_curveFun_ymodel("BLS", "LCSM"), "does not support")
})

test_that("validate_curveFun_ymodel passes for valid combinations", {
  expect_silent(nlpsem:::validate_curveFun_ymodel("QUAD", "LCSM"))
  expect_silent(nlpsem:::validate_curveFun_ymodel("EXP", "LCSM"))
  expect_silent(nlpsem:::validate_curveFun_ymodel("LIN", "LGCM"))
  expect_silent(nlpsem:::validate_curveFun_ymodel("LIN", NULL))
})
