# Bootstrap Confidence Intervals

# # Statistics of Interest --------------------------------------------------
# get_mean <- function(split, ...) {
#   bt_samp <- analysis(split)
#   theta_i <- mean(bt_samp[["mpg"]])
#   return(theta_i)
# }

# get_trimmed_mean <- function(split, trim, ...) {
#   bt_samp <- analysis(split)
#   theta_i <- mean(bt_samp[["mpg"]], trim = 0.3)
#   return(theta_i)
# }
#
# get_median <- function(split, ...){
#   bt_samp <- analysis(split)
#   theta_i <- median(bt_samp[["mpg"]], ...)
#   return(theta_i)
# }
#
# get_diff_median <- function(splits, ...) {
#   boot_sample <- analysis(splits)
#   theta_i <- median(boot_sample$MonthlyIncome[boot_sample$Gender == "Female"]) -
#     median(boot_sample$MonthlyIncome[boot_sample$Gender == "Male"])
#   return(theta_i)
# }
#
#
# # Helper Functions --------------------------------------------------------
# get_sd <- function(split, variable, ...){
#   bt_samp <- analysis(split)
#   theta_sd <- sd(bt_samp[[variable]])
#   return(theta_sd)
# }
#
#
# get_length <- function(split, variable,...){
#   bt_samp <- analysis(split)
#   bt_sample_size <- length(bt_samp[[variable]])
#   return(bt_sample_size)
# }
#
#
# # High-Level API ----------------------------------------------------------
# boot_ci <- function(bt_resamples, statistic, variable, method = "percentile", level = 0.95, ...) {
#   alpha = 1 - level
#   variables = variable
#   apparent_resample <- bt_resamples$splits[[length(bt_resamples$splits)]]
#   apparent_df <- as.data.frame(apparent_resample)
#   apparent_vals <- apparent_df[[variables]]
#
#   # TO-DO
#   # Sure, here I care abut the mean. But how do I generalize for cases
#   # where the mean is NOT theta_obs (the statistic of interest)?
#   #
#   # use invoke() function at the expense of making boot_ci() call
#   # annoying with yet ANOTHER added parameter to write in the call
#   theta_obs <- invoke(statistic, apparent_vals)
#   # theta_obs <- mean(apparent_vals)
#
#   if (method == "percentile") {
#     # don't do detailed computations in this block, just call
#     results <- boot_ci_perc(bt_resamples, alpha, apparent_vals, ...)
#     return(results)
#   } # interesting side-effect..doesn't get called if not possible to calculate?
#   # or is it just a bug?
#   if (method == "pivot-t"){
#     results <- boot_ci_t(bt_resamples, alpha, apparent_vals, theta_obs, ...)
#     return(results)
#   }
#   if (method == "bca"){
#     results <- boot_ci_bca(bt_resamples, alpha, apparent_vals, ...)
#     return(results)
#   }
#   if (method == "abc"){
#     results <- boot_ci_abc(bt_resamples, alpha, apparent_vals, ...)
#     return(results)
#
#   }
# }
#
#
# # Low-Level API -----------------------------------------------------------
# # TO-DO return tibble with upper lower alpha
# boot_ci_perc <- function(bt_resamples, alpha, data) {
#   ci_perc <-  quantile(bt_resamples$theta_i, probs = c(alpha/2, 1-alpha/2))
#   return(ci_perc)
# }
#
# # TO-DO return tibble with upper lower alpha
# boot_ci_t <- function(bt_resamples, alpha, data, theta_obs) {
#   theta_sd <- map_dbl(bt_resamples$splits, get_sd)
#   bt_sample_size <- map_dbl(bt_resamples$splits, get_length)
#   # theta_obs = bt_resamples$theta_obs[1]
#   theta_b <- (bt_resamples$theta_i-theta_obs)/theta_sd/sqrt(bt_sample_size)
#   bootstrap_t <-quantile(theta_b, probs = c(alpha/2, 1-alpha/2))
#   ci_t <-  theta_obs + bootstrap_t * sd(data) / sqrt(length(data))
#   return (ci_t)
# }
#
#
# # TO-DO return tibble with upper lower alpha
# boot_ci_bca <- function(bt_resamples, alpha, data){
#   theta_hat = mean(bt_resamples$theta_i)
#
#   ### Estimating Z0:
#   po = mean(bt_resamples$theta_i <= theta_hat)
#   Z0 = qnorm(po)
#   Za = qnorm(1-alpha/2)
#
#   # // TO-DO I'm sure you can do LOO sampling using rsample()...akin to LOO CV
#   # loo_rsets <- loo_cv(as_tibble(data))
#   # loo_df <- loo_rsets %>%
#   #   mutate(theta_i =  mean())
#   # Error: C stack usage  7970176 is too close to the limit
#
#   leave_one_out_theta = sapply(1:length(data), function(i){
#     leave_out_data = data[-i] # leave out the ith observation
#     theta_i = mean(leave_out_data)
#     return(theta_i)
#   })
#
#   theta_minus_one = mean(leave_one_out_theta)
#   a = sum( (theta_minus_one - leave_one_out_theta)^3)/( 6 *(sum( (theta_minus_one - leave_one_out_theta)^2))^(3/2) )
#
#   Zu = (Z0+Za)/(1-a*(Z0+Za)) + Z0 # upper limit for Z
#   Zl = (Z0-Za)/(1-a*(Z0-Za)) + Z0 # Lower limit for Z
#   lower_percentile = pnorm(Zl,lower.tail = TRUE) # percentile for Z
#   upper_percentile = pnorm(Zu,lower.tail = TRUE) # percentile for Z
#   ci_bca = as.numeric(quantile(bt_resamples$theta_i, c(lower_percentile,upper_percentile))) # putting those percentiles in place of alpha/2, 1-alpha/
#   return(ci_bca)
# }
#
#
# boot_ci_abc <- function(bt_resamples, alpha, data){
#   paste("Approximate Bootstrap CI Results \n a method of approximating the BCa Intervals
#         \n analytically, without using any Monte Carlo replications at all. \n
#         The S function is called: abcnon(x, tt)")
# }


# get_tmean <- function(x)
#   map_dbl(x,
#           function(x)
#             mean(analysis(x)[["Sepal.Width"]], trim = 0.1)
#   )
# get_tmean(bt$splits)


# things to check for:
# missing data
# all same estimates
# make sure that bt_resamples is a boostrap object (inherits)
# make sure theta_obs is not NA
# make sure that z_pntl has two unique values
# check against rand normal data and standard CI
boot_ci_t <- function(bt_resamples, var, alpha, data = NULL, theta_obs) {
  theta_obs <- theta_obs[[var]]
  theta_se <- sd(bt_resamples[[var]], na.rm = TRUE)/
    sqrt(sum(!is.na((bt_resamples[[var]]))))
  z_dist <- (bt_resamples[[var]] - theta_obs)/theta_se
  z_pntl <- quantile(z_dist, probs = c(alpha/2, 1 - (alpha)/2), na.rm = TRUE)
  ci <- theta_obs + z_pntl * theta_se
  tibble(
    lower = ci[1],
    upper = ci[2],
    alpha = alpha,
    method = "bootstrap-t"
  )
}










