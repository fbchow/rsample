#' Bootstrap Confidence Intervals
#'
#' @details
#' Calculate boostrap confidence intervals for a statistic of interest
#'
#'  @importFrom dplyr mutate
#'  @importFrom stats sd
#'  @importFrom dplyr as_tibble
#'  @export
boot_ci_t <- function(bt_resamples, alpha, data = NULL) {

  est <- bt_resamples %>% filter(id == "Apparent") %>% pull(wt_est)
  se <- bt_resamples %>% filter(id == "Apparent") %>% pull(wt_var) %>% sqrt()
  pctl <- bt_resamples %>% filter(id != "Apparent") %>% pull(Z) %>%
  quantile(probs = c(alpha / 2, 1 - alpha / 2)) %>% rev() %>% unname()
  student_ci <- est - pctl * se

  if (all(is.na(est)))
    stop("All bootstrap resample estimates (est) are missing values.", call. = FALSE)

   if (se == 0 | se == Inf)
    stop("Your standard error (se) is 0 or infinity.", call. = FALSE)


   tibble(
    lower = student_ci[1],
    upper = student_ci[2],
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


boot_ci_bca <- function(bt_resamples, func, alpha, data = NULL){

  if (nrow(bt_resamples) < 1000)
    warning("Recommend at least 1000 bootstrap resamples.", call. = FALSE)

  # TODO then write a test case for that
  # if(apparent != TRUE)
  #   warning("Please set apparent = TRUE in bootsraps()")

  #   get_loo_est <-  function(dat) {
  #   lm_fit <- lm(mpg ~ ., data = dat)
  #   coef(lm_fit)["disp"]
  # }
  # theta_hat <- get_loo_est(dat)
  theta_hat <- bt_resamples[["original"]][1]


  ### Estimating Z0 bias-correction
  po <- mean(bt_resamples[["wt_est"]] <= theta_hat, na.rm = TRUE)
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
  ci_bca <- as.numeric(quantile(bt_resamples[["Z"]], c(lower_percentile, upper_percentile))) # use percentiles in place of (alpha / 2) and  (1 - alpha / 2)

  tibble(
  lower = ci_bca[1],
  upper = ci_bca[2],
  alpha = alpha,
  method = "BCa"
  )
}
