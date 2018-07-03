context("boot_ci")
library(rsample)
library(testthat)
library(purrr)
library(tibble)
library(dplyr)


# Try re-writing this example again later --------------------------------
# set.seed(888)
# disp_effect <- function(dat) {
#   lm_fit <- lm(mpg ~ ., data = dat)
#   coef(lm_fit)["disp"]
# }

# bt_splits <- bootstraps(mtcars, times = 5000, apparent = TRUE)
# bt_splits <- bt_splits %>%
#   as_tibble() %>%
#   mutate(
#     model = map(splits, function(x) lm(mpg ~ ., data = analysis(x))),
#     wt_est = map_dbl(model, function(x) coef(x)["wt"]),
#     wt_var = map_dbl(model, function(x) vcov(x)["wt", "wt"]))
# bt_splits <- bt_splits %>%
#   mutate(
#     original = bt_splits %>% filter(id == "Apparent") %>% pull(wt_est),
#     Z = (wt_est - original) / sqrt(wt_var)
#   )

# results_perc <- rsample:::boot_ci_perc(bt_splits,
#                                        stat = "Z",
#                                        alpha = 0.05,
#                                        data = NULL)

# results_t <- rsample:::boot_ci_t(bt_splits,
#                                  stat = "wt_est",
#                                  stat_var = "wt_var",
#                                  alpha = 0.05,
#                                  data = NULL,
#                                  theta_obs = "original",
#                                  var_obs = " ")


# results_bca <- rsample:::boot_ci_bca(bt_splits,
#                                      theta_obs = "original",
#                                      stat = "wt_est",
#                                      func = disp_effect,
#                                      Z = "Z",
#                                      alpha = 0.05,
#                                      data = NULL)




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
      lower = ttest$conf.int[1],
      upper = ttest$conf.int[2],
      alpha = 0.05,
      method = "t-test"
    )

    set.seed(888)
    get_mean <- function(dat){
      mean(dat$rand_nums, na.rm = TRUE)
    }

    bt_norm <- bootstraps(x, times = 1000, apparent = TRUE) %>%
      dplyr::mutate(tmean = map_dbl(splits, function(x) get_mean(analysis(x))))

    bt_norm$original <- mean(x$rand_nums, na.rm=TRUE)


    # results_mean_boot_perc <- rsample:::boot_ci_perc(
    #   bt_norm,
    #   stat = "tmean",
    #   alpha = 0.05,
    #   data = NULL
    # )

    # test pivot-t confidemce interval after refactoring
    #   results_mean_boot_t <- rsample:::boot_ci_t(
    #     bt_norm,
    # )


    results_mean_boot_bca <- rsample:::boot_ci_bca(
      bt_norm,
      theta_obs = "original",
      stat = "tmean",
      func = get_mean,
      Z = "tmean",
      alpha = 0.05,
      data = NULL
    )

    # expect_equal(results_ttest$lower, results_boot_t$lower, tolerance = 0.01)
    # expect_equal(results_ttest$upper, results_boot_t$upper, tolerance = 0.01)
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
})




context("boot_ci: Sufficient Number of Bootstrap Resamples")

get_median <- function(dat){
  median(dat$Sepal.Length, na.rm = TRUE)
}

set.seed(888)
  bt_one <- bootstraps(iris, apparent = TRUE, times = 1) %>%
    dplyr::mutate(median_len = map_dbl(splits, function(x) get_median(analysis(x))))


test_that("At least B=1000 replications needed to sufficiently reduce Monte Carlo sampling Error for BCa method",{
  expect_warning(
    rsample:::boot_ci_bca(
      bt_one,
      theta_obs = "original",
      stat = "median_len",
      func = get_median,
      Z = "median_len",
      alpha = 0.05,
      data = NULL
    )
  )
})
