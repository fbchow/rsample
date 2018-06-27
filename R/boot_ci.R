#' Bootstrap Confidence Intervals
#'
#' @details
#' Calculate boostrap confidence intervals for a statistic of interest
#'
#'  @importFrom dplyr mutate
#'  @importFrom stats sd
#'  @export
boot_ci_t <- function(bt_resamples, stat, alpha, data = NULL, theta_obs) {

  theta_obs <- theta_obs[[stat]]

  if (all(is.na(theta_obs)))
    stop("All statistics (theta_obs) are missing values.", call. = FALSE)

  theta_se <- sd(bt_resamples[[stat]], na.rm = TRUE) / sqrt(sum(!is.na((bt_resamples[[stat]]))))

  if (theta_se == 0 | theta_se == Inf)
    stop("Your standard error (theta_se) is 0 or infinity.", call. = FALSE)

  z_dist <- (bt_resamples[[stat]] - theta_obs) / theta_se
  z_pntl <- quantile(z_dist, probs = c(alpha / 2, 1 - (alpha) / 2), na.rm = TRUE)
  ci <- theta_obs + z_pntl * theta_se

  tibble(
    lower = ci[1],
    upper = ci[2],
    alpha = alpha,
    method = "bootstrap-t"
  )
}

boot_ci_perc <- function(bt_resamples, stat, alpha, data = NULL, theta_obs) {
  z_dist <- bt_resamples[[stat]]

  if (all(is.na(z_dist)))
  stop("All statistics (z_dist) are missing values.", call. = FALSE)

  if (0 < alpha && alpha > 1)
  stop("Your significance level (alpha) is unreasonable.", call. = FALSE)

  ci <- quantile(z_dist, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  tibble(
    lower = ci[1],
    upper = ci[2],
    alpha = alpha,
    method = "percentile"
  )
}


boot_ci_bca <- function(bt_resamples, stat, alpha, data = NULL){

  # TODO then write a test case for that
  if (nrow(bt_resamples) < 1000)
    warning("Recommend at least 1000 bootstrap resamples.", call. = FALSE)

  # TODO then write a test case for that
  # if(apparent != TRUE)
  #   warning("Please set apparent = TRUE in bootsraps()")

  dat <- bt_resamples %>%
    filter(id == "Apparent") %>%
    analysis() %>%
    pluck("splits", 1, "data")
# dat <- dat[["splits"]][[1]][["data"]]

# run this median test again
  # median_difference <- median(dat$MonthlyIncome[dat$Gender == "Female"]) - median(dat$MonthlyIncome[dat$Gender == "Male"])


  # TODO
  # still want this to be generalised
  # want to call in get_theta_i for beta coefficient from original data set
  # theta_hat <- mean(bt_resamples[[stat]], na.rm = TRUE)
  get_theta_i <-  function(dat) {
    lm_fit <- lm(mpg ~ ., data = dat)
    coef(lm_fit)["disp"]
  }
  theta_hat <- get_theta_i(dat)


  ### Estimating Z0 bias-correction
  po <- mean(bt_resamples[[stat]] <= theta_hat)
  Z0 <- qnorm(po)
  Za <- qnorm(1 - alpha / 2)

  # leave_one_out_theta = sapply(1:length(data), function(i){
  #   leave_out_data = data[-i] # leave out the ith observation
  #   theta_i = mean(leave_out_data)
  # return(theta_i)    # returns a vector of means. mean of each bootstrap resample.
  # })

  # TODO
  # generalise to different funtions
  # want to call in the `disp_effect` function here
  leave_one_out_theta <- loo_cv(dat) %>%
    analysis() %>%
    pluck("splits")
    # mutate(theta_i = get_theta_i(splits))

  # `func` is name of bca argument for the orginal function
  leave_one_out_theta <- loo_cv(dat) %>%
    mutate(loo_est = map_dbl(splits, function(x, f) f(analysis(x)), f = func))

  theta_minus_one <- mean(leave_one_out_theta$theta_i)
  a <- sum( (theta_minus_one - leave_one_out_theta$theta_i) ^ 3) / ( 6 * (sum( (theta_minus_one - leave_one_out_theta$theta_i) ^ 2)) ^ (3 / 2) )

  Zu <- (Z0 + Za) / ( 1 - a * (Z0 + Za)) + Z0 # upper limit for Z
  Zl <- (Z0 - Za) / (1 - a * (Z0 - Za)) + Z0 # Lower limit for Z
  lower_percentile <- pnorm(Zl, lower.tail = TRUE) # percentile for Z
  upper_percentile <- pnorm(Zu, lower.tail = TRUE) # percentile for Z
  ci_bca <- as.numeric(quantile(bt_resamples[[stat]], c(lower_percentile, upper_percentile))) # use percentiles in place of (alpha / 2) and  (1 - alpha / 2)

  tibble(
  lower = ci_bca[1],
  upper = ci_bca[2],
  alpha = alpha,
  method = "BCa"
  )
}


