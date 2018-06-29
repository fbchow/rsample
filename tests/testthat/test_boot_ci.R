context("boot_ci")
library(rsample)
library(testthat)
library(purrr)
library(tibble)
library(dplyr)


# Example Code -------------------------------------------------------

# boostrap single statistic k = 1
# get_tmean <- function(x)
#   map_dbl(x,
#           function(x)
#             mean(analysis(x)[["Sepal.Width"]], trim = 0.1))

set.seed(55)
bt_splits <- bootstraps(mtcars, times = 5000, apparent = TRUE) %>%
  as_tibble() %>%
  mutate(
    model = map(splits, function(x) lm(mpg ~ ., data = analysis(x))),
    wt_est = map_dbl(model, function(x) coef(x)["wt"]),
    wt_var = map_dbl(model, function(x) vcov(x)["wt", "wt"]),
    original = bt_splits %>% filter(id == "Apparent") %>% pull(wt_est),
    Z = (wt_est - original)/sqrt(wt_var)
  )



# TODO try diff of medians
data("attrition")
median_diff <- function(splits) {
  x <- analysis(splits)
  median(x$MonthlyIncome[x$Gender == "Female"]) -
    median(x$MonthlyIncome[x$Gender == "Male"])
}

# stat_func <- function(splits) {
#   x <- analysis(splits)
#   median(x$MonthlyIncome[x$Gender == "Female"]) -
#     median(x$MonthlyIncome[x$Gender == "Male"])
# }


set.seed(353)
boot_resamples <- bootstraps(attrition, times = 1000, apparent = TRUE)
boot_resamples$wage_diff <- map_dbl(boot_resamples$splits, median_diff)
boot_resamples

#
# results_median <- rsample:::boot_ci_bca(
#   bt_resamples = boot_resamples %>% dplyr::filter(id != "Apparent"),
#   stat = "wage_diff",
#   stat_func = median_diff,
#   alpha = 0.05,
#   var = "MonthlyIncome",
#   theta_obs = boot_resamples %>% dplyr::filter(id == "Apparent")
# )
#
# results_median



# example to test with: getting a regression coef for one predictor
disp_effect <- function(dat) {
  lm_fit <- lm(mpg ~ ., data = dat)
  coef(lm_fit)["disp"]
}

set.seed(55)
bt_splits <- bootstraps(mtcars, times = 20, apparent = TRUE) %>%
  mutate(beta = map_dbl(splits, function(x) disp_effect(analysis(x))))


bt_splits <- bootstraps(mtcars, times = 5000, apparent = TRUE) %>%
  as_tibble() %>%
  mutate(
    model = map(splits, function(x) lm(mpg ~ ., data = analysis(x))),
    wt_est = map_dbl(model, function(x) coef(x)["wt"]),
    wt_var = map_dbl(model, function(x) vcov(x)["wt", "wt"]))

bt_splits <- bt_splits %>%
  mutate(
    original = bt_splits %>% filter(id == "Apparent") %>% pull(wt_est),
    Z = (wt_est - original)/sqrt(wt_var)
  )

boot_ci_t(bt_splits, alpha = 0.5, data = NULL)






results_coeff <- rsample:::boot_ci_bca(
  bt_resamples = bt_splits,
  stat = "beta",
  alpha = 0.05
)





set.seed(888)
bt_one <- bootstraps(iris, apparent = TRUE, times = 1) %>%
  dplyr::mutate(tmean = get_tmean(splits))

bt <- bootstraps(iris, apparent = TRUE, times = 1000) %>%
  dplyr::mutate(tmean = get_tmean(splits))
#
# bt_lm <- boostraps(iris, apparent = TRUE, times = 1000) %>%
#   dplyr::mutate(tmean = get_tmean(splits))

results_t <- rsample:::boot_ci_t(
  bt_resamples = bt %>% dplyr::filter(id != "Apparent"),
  stat = "tmean",
  alpha = 0.05,
  theta_obs = bt %>% dplyr::filter(id == "Apparent")
)

results_percentile <- rsample:::boot_ci_perc(
  bt_resamples = bt_one %>% dplyr::filter(id != "Apparent"),
  stat = "tmean",
  alpha = 0.05,
  theta_obs = bt_one %>% dplyr::filter(id == "Apparent")
)

results_bca <- rsample:::boot_ci_bca(
  bt_resamples = bt %>% dplyr::filter(id != "Apparent"),
  stat = "tmean",
  alpha = 0.05,
  var = "Sepal.Width",
  theta_obs = bt %>% dplyr::filter(id == "Apparent")
)


context("boot_ci: Sufficient Number of Bootstrap Resamples")
test_that("throw warning if theta_se equals 0 or infinity", {
  set.seed(888)
  bt_one <- bootstraps(iris, apparent = TRUE, times = 1) %>%
    dplyr::mutate(tmean = get_tmean(splits))

  expect_error(
    rsample:::boot_ci_t(
      bt_resamples = bt_one %>% dplyr::filter(id != "Apparent"),
      stat = "tmean",
      alpha = 0.05,
      theta_obs = bt_one %>% dplyr::filter(id == "Apparent")
    )
  )
})
test_that("At least B=1000 replications needed to sufficiently reduce Monte Carlo sampling Error for BCa method",{
  expect_error(
    rsample:::boot_ci_bca(
      bt_resamples = bt_one %>% dplyr::filter(id != "Apparent"),
      stat = "tmean",
      alpha = 0.05,
      theta_obs = bt_one %>% dplyr::filter(id == "Apparent")
    )
  )
})


test_that('z_pntl has two unique values', {
  expect_false(results_t$lower == results_t$upper)
  expect_true(results_percentile$lower == results_percentile$upper)
})

test_that('bootstrap resample estimates are unique',{
  times <- 1
  bt_same <- bootstraps(iris, apparent = TRUE, times = times) %>%
    dplyr::mutate(tmean = rep(3, times + 1))
  expect_error(
    rsample:::boot_ci_t(
      bt_resamples = bt_same %>% dplyr::filter(id != "Apparent"),
      stat = "tmean",
      alpha = 0.05,
      theta_obs = bt_same %>% dplyr::filter(id == "Apparent")
    )
  )
})


context("boot_ci Prompt Errors: Too Many Missing Values")
test_that('upper & lower confidence interval does not contain NA', {
  iris_na<- iris
  iris_na$Sepal.Width[c(1, 51, 101)] <- NA

  set.seed(888)
  bt_na <- bootstraps(iris_na, apparent = TRUE, times = 10000) %>%
    dplyr::mutate(tmean = rep(NA_real_, 10001))

  expect_error(
    rsample:::boot_ci_t(
      bt_resamples = bt_na %>% dplyr::filter(id != "Apparent"),
      stat = "tmean",
      alpha = 0.05,
      theta_obs = bt_na %>% dplyr::filter(id == "Apparent")
    )
  )
  expect_error(
    rsample:::boot_ci_perc(
      bt_resamples = bt_na %>% dplyr::filter(id != "Apparent"),
      stat = "tmean",
      alpha = 0.05,
      theta_obs = bt_na %>% dplyr::filter(id == "Apparent")
    )
  )

})

test_that('theta_obs is not NA', {
  expect_equal(sum(is.na(bt_one$tmean)), 0)
})

test_that('bt_resamples is a bootstrap object', {
  expect_equal(class(bt_one)[1], "bootstraps")
})


test_that('alpha is a reasonable level of significance', {
  expect_error(
    rsample:::boot_ci_perc(
      bt_resamples = bt_na %>% dplyr::filter(id != "Apparent"),
      stat = "tmean",
      alpha = 5,
      theta_obs = bt_na %>% dplyr::filter(id == "Apparent")
    )
  )
})


context("boot_ci Check Against Standard Confidence Interval")
test_that(
  'Bootstrap estimate of mean is close to estimate of mean from normal distribution',
  {
    set.seed(888)
    n <- 10000
    mean <- 10
    sd <- 1

    rand_nums <- rnorm(n, mean, sd)
    x <- as.data.frame(rand_nums)

    ttest <- t.test(x)

    results_ttest <- tibble(
      lower = ttest$conf.int[1],
      upper = ttest$conf.int[2],
      alpha = 0.05,
      method = "bootstrap-t"
    )

    get_mean <- function(y)
      map_dbl(y,
              function(y)
                mean(analysis(y)[["rand_nums"]], na.rm = TRUE))

    bt_norm <- bootstraps(x, times = 1000, apparent = TRUE) %>%
      dplyr::mutate(tmean = get_mean(splits))

    results_boot_t <- rsample:::boot_ci_t(
      bt_resamples = bt_norm %>% dplyr::filter(id != "Apparent"),
      stat = "tmean",
      alpha = 0.05,
      theta_obs = bt_norm %>% dplyr::filter(id == "Apparent")
    )

    results_boot_bca <- rsample:::boot_ci_bca(
      bt_resamples = bt_norm %>% dplyr::filter(id != "Apparent"),
      stat = "tmean",
      alpha = 0.05,
      var = "rand_nums",
      theta_obs = bt_norm %>% dplyr::filter(id == "Apparent")
    )

      expect_equal(results_ttest$lower, results_boot_t$lower, tolerance = 0.01)
      expect_equal(results_ttest$upper, results_boot_t$upper, tolerance = 0.01)
      expect_equal(results_ttest$lower, results_boot_bca$lower, tolerance = 0.01)
      expect_equal(results_ttest$upper, results_boot_bca$upper, tolerance = 0.01)

  }
)



