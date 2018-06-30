#' Bootstrap Confidence Intervals
#'
#' @details
#' Calculate boostrap confidence intervals for a statistic of interest
#'
#'  @importFrom dplyr mutate
#'  @importFrom stats sd
#'  @importFrom dplyr as_tibble
#'  @export
boot_ci_perc <- function(bt_resamples, Z, alpha, data = NULL, theta_obs) {
  z_dist <- bt_resamples[[Z]]

  if (all(is.na(z_dist)))
  stop("All statistics (z_dist) are missing values.", call. = FALSE)

  if (0 < alpha && alpha > 1)
  stop("Your significance level (alpha) is unreasonable.", call. = FALSE)

  ci <- quantile(z_dist, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  tibble(
    lower = min(ci),
    upper = max(ci),
    alpha = alpha,
    method = "percentile"
  )
}
