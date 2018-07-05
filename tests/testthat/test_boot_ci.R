context("boot_ci")
library(rsample)
library(testthat)
library(purrr)
library(tibble)
library(dplyr)

context("boot_ci() Check Against Standard Confidence Interval")
test_that(
  'Bootstrap estimate of mean is close to estimate of mean from normal distribution',
  {
    n <- 10000
    mean <- 10
    sd <- 1

    set.seed(888)
    rand_nums <- rnorm(n, mean, sd)
    x <- as.data.frame(rand_nums)

    ttest <- t.test(x)

    results_ttest <- tibble(
      lower = min(ttest$conf.int),
      upper = max(ttest$conf.int),
      alpha = 0.05,
      method = "t-test"
    )

    set.seed(888)
    get_mean <- function(dat){
      mean(dat$rand_nums, na.rm = TRUE)
    }

    get_sd <- function(dat){
      sd(dat$rand_nums, na.rm = TRUE)
    }

    bt_norm <- bootstraps(x, times = 1000, apparent = TRUE) %>%
      dplyr::mutate(tmean = map_dbl(splits, function(x) get_mean(analysis(x))),
                    tmean_var = map_dbl(splits, function(x) get_sd(analysis(x))))

      results_mean_boot_t <- rsample:::boot_ci_t(
        bt_norm,
        stat = "tmean",
        stat_var = "tmean_var",
        alpha = 0.05,
        data = NULL
    )

    results_mean_boot_bca <- rsample:::boot_ci_bca(
      bt_norm,
      stat = "tmean",
      func = get_mean,
      alpha = 0.05,
      data = NULL
    )

    expect_equal(results_ttest$lower, results_mean_boot_t$lower, tolerance = 0.01)
    expect_equal(results_ttest$upper, results_mean_boot_t$upper, tolerance = 0.01)
    expect_equal(results_ttest$lower, results_mean_boot_bca$lower, tolerance = 0.01)
    expect_equal(results_ttest$upper, results_mean_boot_bca$upper, tolerance = 0.01)
  }
)





context("boot_ci() Prompt Errors: Too Many Missing Values")
test_that('Upper & lower confidence interval does not contain NA', {
  iris_na <- iris
  iris_na$Sepal.Width[c(1, 51, 101)] <- NA

  set.seed(888)
  bt_na <- bootstraps(iris_na, apparent = TRUE, times = 10000) %>%
    dplyr::mutate(tmean = rep(NA_real_, 10001))

  expect_error(
    rsample:::boot_ci_perc(
      bt_na,
      stat = "tmean",
      alpha = 0.05,
      data = NULL
    )
  )

  expect_error(
    rsample:::boot_ci_t(
      bt_na,
      stat = "tmean",
      stat_var = "tmean_var",
      alpha = 0.05,
      data = NULL
      )
  )

  expect_error(
    rsample:::boot_ci_bca(
      bt_na,
      stat = "tmean",
      stat_var = get_mean,
      alpha = 0.05,
      data = NULL
    )
  )

})




context("boot_ci: Sufficient Number of Bootstrap Resamples")

get_median <- function(dat){
  median(dat$Sepal.Length, na.rm = TRUE)
}

set.seed(888)
  bt_one <- bootstraps(iris, apparent = TRUE, times = 1) %>%
    dplyr::mutate(median_sepal = map_dbl(splits, function(x) get_median(analysis(x))))


test_that("At least B=1000 replications needed to sufficiently reduce Monte Carlo sampling Error for BCa method",{
  expect_warning(
    rsample:::boot_ci_bca(
      bt_one,
      stat = "median_sepal",
      func = get_median,
      alpha = 0.05,
      data = NULL
    )
  )
})


test_that('theta_obs is not NA', {
  expect_equal(sum(is.na(bt_one$median_sepal)), 0)
})

test_that('bt_resamples is a bootstrap object', {
  expect_equal(class(bt_one)[1], "bootstraps")
})
