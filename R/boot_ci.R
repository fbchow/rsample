#' Bootstrap Confidence Intervals
#'
#' @details
#' Calculate boostrap confidence intervals for a statistic of interest
#'
#'  @export

#' @importFrom dplyr mutate
#' @importFrom stats sd quantile pnorm
#' @importFrom dplyr as_tibble
#' @export
boot_ci_t <- function(bt_resamples, stat, stat_var, alpha = 0.05, data = NULL, theta_obs, var_obs) {

  theta_obs <- theta_obs[[stat]]
  var_obs <- var_obs[[stat_var]]

  if (all(is.na(theta_obs)))
    stop("All statistics (theta_obs) are missing values.", call. = FALSE)

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



#' @export
boot_ci_perc <- function(bt_resamples, stat, alpha = 0.05, data = NULL, theta_obs = NULL) {

  if (all(is.na(bt_resamples[[stat]])))
    stop("All statistics (", stat, ") are missing values.", call. = FALSE)

  if (0 < alpha && alpha > 1)
    stop("Your significance level (alpha) is unreasonable.", call. = FALSE)

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

#' @export
boot_ci_bca <- function(bt_resamples, theta_obs, stat, func, Z, alpha = 0.05, data = NULL){

  if (nrow(bt_resamples) < 1000)
    warning("Recommend at least 1000 bootstrap resamples.", call. = FALSE)

  theta_hat <- bt_resamples[[theta_obs]][1]


  ### Estimating Z0 bias-correction
  po <- mean(bt_resamples[[stat]] <= theta_hat, na.rm = TRUE)
  Z0 <- qnorm(po)
  Za <- qnorm(1 - alpha / 2)

  # `func` is name of bca argument for the orginal function
  leave_one_out_theta <- loo_cv(bt_resamples %>%
                                  filter(id == "Apparent") %>%
                                  pluck("splits", 1, "data")) %>%
    mutate(loo_est = map_dbl(splits, function(x, f) f(analysis(x)), f = func))

  theta_minus_one <- mean(leave_one_out_theta$loo_est, na.rm = TRUE)
  a <- sum( (theta_minus_one - leave_one_out_theta$loo_est) ^ 3) / ( 6 * (sum( (theta_minus_one - leave_one_out_theta$loo_est) ^ 2)) ^ (3 / 2) )

  Zu <- (Z0 + Za) / ( 1 - a * (Z0 + Za)) + Z0 # upper limit for Z
  Zl <- (Z0 - Za) / (1 - a * (Z0 - Za)) + Z0 # Lower limit for Z
  lower_percentile <- pnorm(Zl, lower.tail = TRUE) # percentile for Z
  upper_percentile <- pnorm(Zu, lower.tail = TRUE) # percentile for Z
  ci_bca <- as.numeric(quantile(bt_resamples[[Z]], c(lower_percentile, upper_percentile))) # use percentiles in place of (alpha / 2) and  (1 - alpha / 2)

  tibble(
  lower = min(ci_bca),
  upper = max(ci_bca),
  alpha = alpha,
  method = "BCa"
  )
}
