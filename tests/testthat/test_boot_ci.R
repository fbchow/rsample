context("boot_ci")
library(rsample)
library(testthat)
library(purrr)
library(tibble)
library(dplyr)

# EX 1: estimating regression coeff for one predictor
disp_effect <- function(dat) {
  lm_fit <- lm(mpg ~ ., data = dat)
  coef(lm_fit)["disp"]
}

# set.seed(55)
bt_splits <- bootstraps(mtcars, times = 5000, apparent = TRUE)
bt_splits <- bt_splits %>%
  as_tibble() %>%
  mutate(
    model = map(splits, function(x) lm(mpg ~ ., data = analysis(x))),
    wt_est = map_dbl(model, function(x) coef(x)["wt"]),
    wt_var = map_dbl(model, function(x) vcov(x)["wt", "wt"]))
bt_splits <- bt_splits %>%
  mutate(
    original = bt_splits %>% filter(id == "Apparent") %>% pull(wt_est),
    Z = (wt_est - original) / sqrt(wt_var)
  )

results_perc <- rsample:::boot_ci_perc(bt_splits, Z = "Z", alpha = 0.05, data = NULL)
