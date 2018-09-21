#' Nonparametric Bootstrap Confidence Intervals
#' @description
#' Calculate bootstrap confidence intervals for a statistic of interest.
#' @param resamples_object An `rsplit` object created by the `bootstraps` function
#' @param ... confidence interval arguments
#' @param method "percentile", "bca", or "student-t"
#' @param alpha level of significance
#' @importFrom rlang eval_tidy quos quo
#' @importFrom tidyselect vars_select
#' @importFrom purrr map_dbl map2_dfr map_dfr
#' @importFrom future plan
#' @importFrom furrr future_map_dfr future_map2_dfr future_map_dbl
#' @export
boot_ci <- function(resamples_object,
                    ...,
                    method = c("percentile", "bca", "student-t"),
                    alpha = 0.05)  {

#TODO - is this right? 
plan(multiprocess)

  method <- match.arg(method)

  other_options <- quos(...)
  other_names <- names(other_options)

  # arguments named at birth of `boot_ci()` function call
  args <- other_options[other_names != ""]

  if (any(names(args) == "func")) {
    bca_names <- names(args)[!(names(args) %in% c("", "stat_var"))]
    bca_args <- args[bca_names]

  } else{
    bca_args <- NULL
  }


  if (any(names(args) == "stat_var")) {
    t_args <- args["stat_var"]

  } else{
    t_args <- NULL
  }


  # arugments with no names at birth `boot_ci()` function call
  predictor_vars <- other_options[other_names == ""]

  # clean, get values, unpack, prepare for evaluation
  predictors <-
    unname(vars_select(names(resamples_object), !!!predictor_vars))


  if (method == "percentile") {
    # We should assume that there are multiple columns being analyzed so we
    # would need to iterate over their names(in `predictors`).
    # Instead of a `for` loop, let's use `map`
    perc_wrapper <- function(stat, bt_resamples, alpha) {
      boot_ci_perc(bt_resamples = bt_resamples,
                   alpha = alpha,
                   stat = stat)
    }


    results <-
      future_map_dfr(predictors,
              perc_wrapper,
              bt_resamples = resamples_object,
              alpha = alpha) %>%
      mutate(variable = predictors)
    return(results)
  }


  if (method == "student-t") {

    var_predictors <-
      unname(vars_select(names(resamples_object), !!!t_args))


    t_wrapper <- function(stat, stat_var, bt_resamples, alpha) {
      boot_ci_t(
        bt_resamples = bt_resamples,
        alpha = alpha,
        stat = stat,
        stat_var = stat_var
      )
    }

    map_expr <-
      future_map2_dfr(predictors,
           var_predictors,
           t_wrapper,
           resamples_object,
           alpha = alpha)

    results <- eval_tidy(map_expr) %>%
      mutate(variable = predictors)
    return(results)

  }


  if (method == "bca") {
    bca_wrapper <- function(stat, bt_resamples, alpha, func, ...) {
      res <- boot_ci_bca(
        bt_resamples = bt_resamples,
        alpha = alpha,
        stat = stat,
        func = func,
        ...
      )
      return(res)
    }

    map_expr <-
      quo(future_map_dfr(
        predictors,
        bca_wrapper,
        resamples_object,
        alpha = alpha,!!!bca_args
      ))
    results <- eval_tidy(map_expr) %>%
      mutate(variable = predictors)
    return(results)

  }
}



#' Student-T Bootstrap Confidence Intervals
#' @description
#' Calculate bootstrap confidence intervals for a statistic of interest.
#' @param bt_resamples An `rsplit` object created by the `bootstraps` function
#' @param stat A statistic of interest
#' @param stat_var The variance of each bootstrap resample
#' @param alpha level of significance
#' @param ... confidence interval arguments
#' @importFrom dplyr mutate
#' @importFrom stats sd quantile pnorm
#' @importFrom dplyr as_tibble
#' @importFrom dplyr last
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom purrr pluck
#' @export
boot_ci_t <- function(bt_resamples, alpha = 0.05, stat, stat_var, ...) {

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
#' @param ... confidence interval arguments
#' @export
boot_ci_perc <- function(bt_resamples, alpha = 0.05, stat, ...) {

  if (all(is.na(bt_resamples[[stat]])))
    stop("All statistics (", stat, ") are missing values.", call. = FALSE)

  # if (0 < alpha && alpha > 1)
    # stop("Your significance level (alpha) is unreasonable.", call. = FALSE)

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
#' @param ... Optional extra arguments to pass to `func`
#' @importFrom furrr future_map_dbl
#' @importFrom furrr future_map_dfr
#' @export
boot_ci_bca <- function(bt_resamples, alpha = 0.05, stat, func, ...){

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

  loo_rs <-
    loo_cv(
      bt_resamples %>%
        dplyr::filter(id == "Apparent") %>%
        pluck("splits", 1, "data")
      )

  # We can't be sure what we will get back from the analysis function.
  # To test, we run on the first LOO data set and see if it is a vector or
  # df
  loo_test <- func(analysis(loo_rs$splits[[1]]), ...)

  if (is.vector(loo_test)) {
    if (length(loo_test) > 1)
      stop("The function should return a single value or a data frame/",
           "tibble.", call. = FALSE)
    leave_one_out_theta <-
      future_map_dbl(loo_rs$splits, function(x) func(analysis(x), ...))
  } else {
    if (!is.data.frame(loo_test))
      stop("The function should return a single value or a data frame/",
           "tibble.", call. = FALSE)
    leave_one_out_theta <-
      future_map_dfr(loo_rs$splits, function(x) func(analysis(x), ...))[[stat]]
  }

  theta_minus_one <- mean(leave_one_out_theta, na.rm = TRUE)
  a <- sum( (theta_minus_one - leave_one_out_theta) ^ 3) /
    ( 6 * (sum( (theta_minus_one - leave_one_out_theta) ^ 2)) ^ (3 / 2) )

  Zu <- (Z0 + Za) / ( 1 - a * (Z0 + Za)) + Z0 # upper limit for Z
  Zl <- (Z0 - Za) / (1 - a * (Z0 - Za)) + Z0 # lower limit for Z
  lower_percentile <- pnorm(Zl, lower.tail = TRUE) # percentile for Z
  upper_percentile <- pnorm(Zu, lower.tail = TRUE) # percentile for Z
  ci_bca <- as.numeric(quantile(bt_resamples[[stat]], c(lower_percentile, upper_percentile)))


  tibble(
  lower = min(ci_bca),
  upper = max(ci_bca),
  alpha = alpha,
  method = "BCa"
  )
}


