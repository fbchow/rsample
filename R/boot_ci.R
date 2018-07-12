
#' Student-T Bootstrap Confidence Intervals
#' @description
#' Calculate bootstrap confidence intervals for a statistic of interest.
#' @param bt_resamples An `rsplit` object created by the `bootstraps` function
#' @param stat A statistic of interest
#' @param stat_var The variance of each bootstrap resample
#' @param alpha
#' @importFrom dplyr mutate
#' @importFrom stats sd quantile pnorm
#' @importFrom dplyr as_tibble
#' @importFrom dplyr last
#' @export
boot_ci_t <- function(bt_resamples, stat, stat_var, alpha = 0.05, data = NULL) {

  theta_obs <- bt_resamples %>% filter(id == "Apparent") %>% pull(!!stat)
  var_obs <- bt_resamples %>% filter(id == "Apparent") %>% pull(!!stat_var)

  if (all(is.na(stat)))
    stop("All statistics (", stat, ") are missing values.", call. = FALSE)

  if(all(is.na(stat_var)))
    stop("All statistics (", stat_var, ") are missing values.", call. = FALSE)

  if (!all(alpha < 1) || !all(alpha > 0))
    stop("All elements of alpha must be in (0,1)")

  if (missing(stat))
    stop("Please specify statistic of interest (stat).")

  z_dist <- (bt_resamples[[stat]] - theta_obs) / sqrt(bt_resamples[[stat_var]])

  z_pntl <- quantile(z_dist, probs = c(alpha / 2, 1 - (alpha) / 2), na.rm = TRUE)
  ci <- theta_obs - z_pntl * sqrt(var_obs)


  tibble(
    lower = min(ci),
    upper = max(ci),
    alpha = alpha,
    method = "bootstrap-t"
  )
}

#' Percentile Method
#' @description
#' Calculate bootstrap confidence intervals for a statistic of interest.
#' @param bt_resamples An `rsplit` object created by the `bootstraps` function
#' @param stat A statistic of interest
#' @param alpha
#' @export
boot_ci_perc <- function(bt_resamples, stat, alpha = 0.05, data = NULL) {

  if (all(is.na(bt_resamples[[stat]])))
    stop("All statistics (", stat, ") are missing values.", call. = FALSE)

  if (0 < alpha && alpha > 1)
    stop("Your significance level (alpha) is unreasonable.", call. = FALSE)

  if (missing(stat))
    stop("Please specify statistic of interest (stat).")

  ci <-
    bt_resamples %>%
    filter(id != "Apparent") %>%
    pull(!!stat) %>%
    quantile(probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)

  tibble(
    lower = min(ci),
    upper = max(ci),
    alpha = alpha,
    method = "percentile"
  )
}

#' BCa Method
#' @description
#' Calculate bootstrap confidence intervals for a statistic of interest.
#' @param bt_resamples An `rsplit` object created by the `bootstraps` function
#' @param stat A statistic of interest
#' @param func A function which when applied to data returns a vector containing the statistic(s) of interest.
#' @param alpha
#' @export
boot_ci_bca <- function(bt_resamples, stat, func, alpha = 0.05, data = NULL){

  if (nrow(bt_resamples) < 1000)
    warning("Recommend at least 1000 bootstrap resamples.", call. = FALSE)

  if (!all(alpha < 1) || !all(alpha > 0))
    stop("All elements of alpha must be in (0,1)")

  if (class(func) != "function")
    stop("Please enter a function to calculate a statistic of interest.")

  if (missing(stat))
    stop("Please specify statistic of interest (stat).")

  if(bt_resamples %>% pull("id") %>% dplyr::last() != "Apparent")
    stop("Please set apparent=TRUE in bootstraps() function")


  theta_hat <- bt_resamples %>% filter(id == "Apparent") %>% pull(!!stat)


  ### Estimating Z0 bias-correction
  # po <- mean(bt_resamples[[stat]] <= theta_hat, na.rm = TRUE)
  po <- mean(bt_resamples %>% pull(!!stat) <= theta_hat, na.rm = TRUE)
  Z0 <- qnorm(po)
  Za <- qnorm(1 - alpha / 2)

  leave_one_out_theta <- loo_cv(bt_resamples %>%
                                  filter(id == "Apparent") %>%
                                  pluck("splits", 1, "data")) %>%
    mutate(loo_est = map_dbl(splits, function(x, f) f(analysis(x)), f = func))

  theta_minus_one <- mean(leave_one_out_theta$loo_est, na.rm = TRUE)
  a <- sum( (theta_minus_one - leave_one_out_theta$loo_est) ^ 3) / ( 6 * (sum( (theta_minus_one - leave_one_out_theta$loo_est) ^ 2)) ^ (3 / 2) )

  Zu <- (Z0 + Za) / ( 1 - a * (Z0 + Za)) + Z0 # upper limit for Z
  Zl <- (Z0 - Za) / (1 - a * (Z0 - Za)) + Z0 # lower limit for Z
  lower_percentile <- pnorm(Zl, lower.tail = TRUE) # percentile for Z
  upper_percentile <- pnorm(Zu, lower.tail = TRUE) # percentile for Z
  # use percentiles in place of (alpha / 2) and  (1 - alpha / 2)
  ci_bca <- as.numeric(quantile(bt_resamples[[stat]], c(lower_percentile, upper_percentile)))


  tibble(
  lower = min(ci_bca),
  upper = max(ci_bca),
  alpha = alpha,
  method = "BCa"
  )
}


