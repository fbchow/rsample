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


context("boot_ci Check Against Standard Confidence Interval")
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
      dplyr::mutate(
        tmean = map_dbl(splits, function(x) get_mean(analysis(x))))

    bt_norm$original <- mean(x$rand_nums, na.rm=TRUE)


    results_mean_boot_perc <- rsample:::boot_ci_perc(
      bt_norm,
      stat = "tmean",
      alpha = 0.05,
      data = NULL
    )

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



