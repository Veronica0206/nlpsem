test_that("validate_columns catches missing y_var columns", {
  dat <- data.frame(T1 = 1, T2 = 2, Y1 = 10, Y2 = 20)
  # Y3 doesn't exist

  expect_error(
    nlpsem:::validate_columns(dat, t_var = "T", y_var = "Y", records = 1:3),
    "Column\\(s\\) not found in data:.*T3.*Y3"
  )
})

test_that("validate_columns catches missing t_var columns", {
  dat <- data.frame(T1 = 1, T2 = 2, Y1 = 10, Y2 = 20, Y3 = 30)
  expect_error(
    nlpsem:::validate_columns(dat, t_var = "T", y_var = "Y", records = 1:3),
    "Column\\(s\\) not found in data:.*T3"
  )
})

test_that("validate_columns passes when all columns exist", {
  dat <- data.frame(T1 = 1, T2 = 2, T3 = 3, Y1 = 10, Y2 = 20, Y3 = 30)
  expect_silent(nlpsem:::validate_columns(dat, t_var = "T", y_var = "Y", records = 1:3))
})

test_that("validate_columns checks growth_TIC columns", {
  dat <- data.frame(T1 = 1, Y1 = 10, ex1 = 0.5)
  expect_error(
    nlpsem:::validate_columns(dat, t_var = "T", y_var = "Y", records = 1, growth_TIC = c("ex1", "ex2")),
    "Column\\(s\\) not found in data:.*ex2"
  )
})

test_that("validate_columns checks TVC columns", {
  dat <- data.frame(T1 = 1, T2 = 2, Y1 = 10, Y2 = 20)
  expect_error(
    nlpsem:::validate_columns(dat, t_var = "T", y_var = "Y", records = 1:2, TVC = "Rs"),
    "Column\\(s\\) not found in data:.*Rs1.*Rs2"
  )
})

test_that("validate_columns checks grp_var column", {
  dat <- data.frame(T1 = 1, Y1 = 10)
  expect_error(
    nlpsem:::validate_columns(dat, t_var = "T", y_var = "Y", records = 1, grp_var = "SEX"),
    "Column\\(s\\) not found in data:.*SEX"
  )
})

test_that("validate_columns works with list records (MGM/MED)", {
  dat <- data.frame(T1 = 1, T2 = 2, R1 = 10, R2 = 20, M1 = 5, M2 = 6)
  expect_silent(
    nlpsem:::validate_columns(dat, t_var = c("T", "T"), y_var = c("R", "M"),
                     records = list(1:2, 1:2))
  )
  expect_error(
    nlpsem:::validate_columns(dat, t_var = c("T", "T"), y_var = c("R", "M"),
                     records = list(1:2, 1:3)),
    "Column\\(s\\) not found in data:.*M3"
  )
})

test_that("validate_columns handles baseline x_type", {
  dat <- data.frame(T1 = 1, Y1 = 10, M1 = 5, age = 25)
  expect_silent(
    nlpsem:::validate_columns(dat, t_var = "T", y_var = "Y", records = list(1, 1),
                     m_var = "M", x_var = "age", x_type = "baseline")
  )
  expect_error(
    nlpsem:::validate_columns(dat, t_var = "T", y_var = "Y", records = list(1, 1),
                     m_var = "M", x_var = "missing_var", x_type = "baseline"),
    "Column\\(s\\) not found in data:.*missing_var"
  )
})

test_that("validate_columns checks cluster_TIC columns", {
  dat <- data.frame(T1 = 1, Y1 = 10, gx1 = 0.5)
  expect_error(
    nlpsem:::validate_columns(dat, t_var = "T", y_var = "Y", records = 1,
                     cluster_TIC = c("gx1", "gx2")),
    "Column\\(s\\) not found in data:.*gx2"
  )
})
