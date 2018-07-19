
#' Student-T Bootstrap Confidence Intervals
#' @description
#' Calculate bootstrap confidence intervals for a statistic of interest.
#' @param bt_resamples An `rsplit` object created by the `bootstraps` function
#' @param stat A statistic of interest
#' @param stat_var The variance of each bootstrap resample
#' @param alpha level of significance
#' @param data NULL
#' @importFrom dplyr mutate
#' @importFrom stats sd quantile pnorm
#' @importFrom dplyr as_tibble
#' @importFrom dplyr last
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom purrr pluck
#' @export
boot_ci_t <- function(bt_resamples, stat, stat_var, alpha = 0.05, data = NULL) {

  theta_obs <- bt_resamples %>%
    dplyr::filter(id == "Apparent") %>%
    dplyr::pull(!!stat)
  var_obs <- bt_resamples %>%
    dplyr::filter(id == "Apparent") %>%
    dplyr::pull(!!stat_var)

  # TODO check if codecovers for this test...
  if (all(is.na(stat)))
    stop("All statistics (", stat, ") are missing values.", call. = FALSE)

  if(all(is.na(stat_var)))
    stop("All statistics (", stat_var, ") are missing values.", call. = FALSE)

  if (!all(alpha < 1) || !all(alpha > 0))
    stop("All elements of alpha must be in (0,1)", call. = FALSE)

  if (missing(stat))
    stop("Please specify statistic of interest (stat).", call. = FALSE)

  # T.b = (theta.i-theta.obs)/(sd(boot.sample)/sqrt(length(boot.sample)))
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
#' @param alpha level of significance
#' @param data NULL
#' @export
boot_ci_perc <- function(bt_resamples, stat, alpha = 0.05, data = NULL) {

  if (all(is.na(bt_resamples[[stat]])))
    stop("All statistics (", stat, ") are missing values.", call. = FALSE)

  if (0 < alpha && alpha > 1)
    stop("Your significance level (alpha) is unreasonable.", call. = FALSE)

  if (missing(stat))
    stop("Please specify statistic of interest (stat).", call. = FALSE)


  ci <-
    bt_resamples %>%
    dplyr::filter(id != "Apparent") %>%
    dplyr::pull(!!stat) %>%
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
#' @param func A function which when applied to data returns a vector containing the statistics of interest.
#' @param alpha level of significance
#' @param data NULL
#' @param ... Optional extra arguments to pass to `func`.
#' @export
boot_ci_bca <- function(bt_resamples, stat, func, alpha = 0.05, data = NULL, ...){

  if (nrow(bt_resamples) < 1000)
    warning("Recommend at least 1000 bootstrap resamples.", call. = FALSE)

  if (!all(alpha < 1) || !all(alpha > 0))
    stop("All elements of alpha must be in (0,1)", call. = FALSE)

  #TODO add call. = false
  if (class(func) != "function")
    stop("Please enter a function to calculate a statistic of interest.", call. = FALSE)

  if (class(bt_resamples)[1] != "bootstraps")
    stop("Please enter a bootstraps sample using the rsample package.", call. = FALSE)

# TODO
  # if (type(func(args)) != data.frame)
  #   stop("Your function", func, "needs to accept a data.frame or tibble as arguments.")

  if (missing(stat))
    stop("Please specify statistic of interest (stat).", call. = FALSE)

  if(bt_resamples %>% pull("id") %>% dplyr::last() != "Apparent")
    stop("Please set apparent=TRUE in bootstraps() function", call. = FALSE)


  theta_hat <- bt_resamples %>%
    dplyr::filter(id == "Apparent") %>%
    dplyr::pull(!!stat)

  ### Estimating Z0 bias-correction
  po <- mean(bt_resamples %>% dplyr::pull(!!stat) <= theta_hat, na.rm = TRUE)
  Z0 <- qnorm(po)
  Za <- qnorm(1 - alpha / 2)

  leave_one_out_theta <- loo_cv(bt_resamples %>%
                                  dplyr::filter(id == "Apparent") %>%
                                  pluck("splits", 1, "data")) %>%
            dplyr::mutate(loo_est =
                            map_dbl(splits, function(x) func(analysis(x), ...)))

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


