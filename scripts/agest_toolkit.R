# ==============================================================================
# AGEST - Age Estimation Toolkit
# ==============================================================================
# R translation of Stata code from Navitainuck et al., 2024
# "Osteological age-at-death estimation in an archaeological sample avoiding 
# age-mimicry"
#
# This toolkit provides functions for estimating age-at-death distributions
# from skeletal remains while correcting for age mimicry bias
# ==============================================================================

# Required libraries
library(tidyverse)
library(splines)

# ==============================================================================
# AUXILIARY FUNCTIONS FOR COMPUTING f(y|a) - PROBABILITY OF OBSERVING 
# CHARACTERISTIC y GIVEN TRUE AGE a
# ==============================================================================

#' Ordinal Regression Probability Model
#' 
#' Computes f(y|a) for ordinal regression models where skeletal characteristics
#' are treated as ordered categories (e.g., stages 1-6 of pubic symphysis aging)
#'
#' @param vary Vector of observed characteristic values (ordinal stages)
#' @param vara Vector of ages
#' @param distf Distribution function ("normal" or "logistic")
#' @param from Minimum characteristic value (e.g., 1)
#' @param to Maximum characteristic value (e.g., 6)
#' @param kappa Vector of threshold parameters (cutpoints between stages)
#' @param trans Optional transformation of age (e.g., "log" or "sqrt")
#' @param sigma Standard deviation parameter (sigma parameterization)
#' @param beta Slope parameter (beta parameterization)
#' @return Vector of probabilities f(y|a)
oregressprob <- function(vary, vara, distf, from, to, kappa, 
                         trans = NULL, sigma = NULL, beta = NULL) {
  
  # Check that either sigma or beta is specified (but not both)
  if (is.null(sigma) && is.null(beta)) {
    stop("Either sigma or beta must be specified")
  }
  
  # Determine which parameterization to use
  if (!is.null(sigma)) {
    param_type <- "sigma"
  } else {
    param_type <- "beta"
  }
  
  # Initialize result vector
  result <- numeric(length(vary))
  
  # Add extreme thresholds at boundaries
  kappa_extended <- c(-1000, kappa, 1000)
  
  # Apply age transformation if specified
  if (!is.null(trans)) {
    # Replace "#" placeholder with actual age variable
    trans_expr <- gsub("#", "vara", trans)
    transa <- eval(parse(text = trans_expr))
  } else {
    transa <- vara
  }
  
  # Select distribution function
  if (distf == "normal") {
    dist_func <- pnorm
  } else if (distf == "logistic") {
    dist_func <- plogis
  } else {
    stop("distf must be 'normal' or 'logistic'")
  }
  
  # For each possible characteristic value, compute probability
  for (y in from:to) {
    # Get upper and lower thresholds for this stage
    idx <- y - from + 2  # Index in extended kappa vector
    upper_kappa <- kappa_extended[idx]
    lower_kappa <- kappa_extended[idx - 1]
    
    # Compute cumulative probabilities based on parameterization
    if (param_type == "sigma") {
      # P(Y <= y | a) - P(Y <= y-1 | a)
      upper_prob <- dist_func((upper_kappa - transa) / sigma)
      lower_prob <- dist_func((lower_kappa - transa) / sigma)
    } else {  # beta parameterization
      upper_prob <- dist_func(upper_kappa - beta * transa)
      lower_prob <- dist_func(lower_kappa - beta * transa)
    }
    
    # Probability of being in this category
    result[vary == y] <- upper_prob[vary == y] - lower_prob[vary == y]
  }
  
  return(result)
}


#' Normal Stratified Probability Model
#' 
#' Computes f(y|a) when y distribution is derived from f(a|y) and P(y)
#' assuming normal distribution of age given characteristic value.
#' This uses Bayes' theorem: f(y|a) ∝ f(a|y) × P(y)
#'
#' @param vary Vector of observed characteristic values
#' @param vara Vector of ages
#' @param from Minimum characteristic value
#' @param to Maximum characteristic value
#' @param mu Vector of mean ages for each characteristic value
#' @param sigma Vector of standard deviations for each characteristic value
#' @param pi Vector of prior probabilities P(y) for each characteristic value
#' @return Vector of probabilities f(y|a)
normalstratprob <- function(vary, vara, from, to, mu, sigma, pi) {
  
  # Initialize result vector
  result <- numeric(length(vary))
  
  # Convert parameter lists to named vectors for easier indexing
  mu_vec <- setNames(mu, from:to)
  sigma_vec <- setNames(sigma, from:to)
  pi_vec <- setNames(pi, from:to)
  
  # For each observation, compute f(a|y) × P(y) for all y values
  sum_probs <- numeric(length(vara))
  
  for (y in from:to) {
    # Likelihood: f(a|y) assuming a ~ Normal(mu[y], sigma[y])
    likelihood <- dnorm(vara, mean = mu_vec[as.character(y)], 
                       sd = sigma_vec[as.character(y)])
    
    # Multiply by prior P(y)
    prob_y <- likelihood * pi_vec[as.character(y)]
    
    # Add to denominator (sum across all y)
    sum_probs <- sum_probs + prob_y
    
    # Store numerator for observations with this y value
    result[vary == y] <- prob_y[vary == y]
  }
  
  # Normalize by dividing by sum (Bayes' theorem)
  result <- result / sum_probs
  
  return(result)
}


#' Regression Density Model
#' 
#' Computes f(y|a) for continuous characteristics following a linear regression
#' model: y = alpha + beta * a + epsilon, where epsilon ~ Normal(0, sigma_eps)
#'
#' @param vary Vector of observed characteristic values (continuous)
#' @param vara Vector of ages
#' @param alpha Intercept parameter
#' @param beta Slope parameter
#' @param sigmaeps Standard deviation of residuals
#' @return Vector of densities f(y|a)
regressdens <- function(vary, vara, alpha, beta, sigmaeps) {
  # The characteristic follows: y ~ Normal(alpha + beta*a, sigma_eps)
  # So density is normal with mean = alpha + beta*a and sd = sigma_eps
  result <- dnorm(vary, mean = alpha + beta * vara, sd = sigmaeps)
  return(result)
}


# ==============================================================================
# VISUALIZATION FUNCTIONS FOR f(y|a) MODELS
# ==============================================================================

#' Visualize Descriptor Model for Ordinal Data
#' 
#' Creates visualization showing how probability of observing each 
#' characteristic stage varies with age
#'
#' @param model_func Function name (e.g., oregressprob)
#' @param params Named list of model parameters
#' @param gridage Vector of ages to evaluate
#' @param ytitle Y-axis title
#' @return ggplot object
showdescriptord <- function(model_func, params, gridage, 
                            ytitle = "Probability") {
  
  # Extract from and to parameters
  from <- params$from
  to <- params$to
  
  # Create grid of all combinations of y values and ages
  data <- expand.grid(
    y = from:to,
    age = gridage
  )
  
  # Remove model-specific parameters before passing to function
  func_params <- params
  func_params$from <- NULL
  func_params$to <- NULL
  
  # Compute density for each combination using the specified model
  data$dens <- do.call(model_func, c(
    list(vary = data$y, vara = data$age, from = from, to = to),
    func_params
  ))
  
  # Compute cumulative probabilities for area plot
  data <- data %>%
    arrange(age, y) %>%
    group_by(age) %>%
    mutate(
      cum = cumsum(dens),
      low = lag(cum, default = 0)
    ) %>%
    ungroup()
  
  # Create stacked area plot showing probability of each stage at each age
  p <- ggplot(data, aes(x = age)) +
    geom_ribbon(aes(ymin = low, ymax = cum, fill = factor(y))) +
    labs(y = ytitle, x = "Age", fill = "Stage") +
    theme_minimal()
  
  return(p)
}


#' Visualize Descriptor Model for Continuous Data
#' 
#' Shows density curves for continuous characteristics at different ages
#'
#' @param model_func Function name
#' @param params Named list of model parameters
#' @param gridchar Vector of characteristic values to evaluate
#' @param gridage Vector of ages to evaluate
#' @return ggplot object
showdescriptcont <- function(model_func, params, gridchar, gridage) {
  
  # Create grid of characteristic values and ages
  data <- expand.grid(
    score = gridchar,
    age = gridage
  )
  
  # Compute density for each combination
  data$dens <- do.call(model_func, c(
    list(vary = data$score, vara = data$age),
    params
  ))
  
  # Create faceted line plots (one panel per age)
  p <- ggplot(data, aes(x = dens, y = score)) +
    geom_line() +
    facet_wrap(~age, nrow = 1) +
    labs(x = "", y = "Score") +
    theme_minimal() +
    theme(axis.text.x = element_blank())
  
  return(p)
}


#' Visualize Model and Observed Data for Ordinal Variable
#' 
#' Side-by-side comparison of model predictions and observed data distribution
#'
#' @param var_data Data frame with the characteristic variable
#' @param var_name Name of the characteristic variable
#' @param model_spec Model specification (function and parameters)
#' @param gridage Ages to evaluate model at
#' @return Combined ggplot object
showdescriptandvarord <- function(var_data, var_name, model_spec, gridage) {
  
  # Extract variable
  var_vec <- var_data[[var_name]]
  
  # Create model visualization
  p1 <- showdescriptord(
    model_func = model_spec$func,
    params = model_spec$params,
    gridage = gridage,
    ytitle = "Probability"
  )
  
  # Create observed data bar plot
  obs_data <- var_data %>%
    filter(!is.na(.data[[var_name]])) %>%
    count(.data[[var_name]]) %>%
    mutate(prob = n / sum(n))
  
  p2 <- ggplot(obs_data, aes(x = factor(.data[[var_name]]), y = prob)) +
    geom_col() +
    labs(x = var_name, y = "Probability", 
         title = paste0(var_name, " (N=", sum(!is.na(var_vec)), ")")) +
    theme_minimal()
  
  # Combine plots
  library(patchwork)
  combined <- p1 + p2
  
  return(combined)
}


# ==============================================================================
# PRIOR DISTRIBUTIONS FOR AGE f(a)
# ==============================================================================

#' Normal Distribution Density
#' 
#' Computes probability density for a normal distribution
#' Used as prior for age distribution f(a)
#'
#' @param age Vector of ages
#' @param mu Mean parameter
#' @param lnsigma Log of standard deviation
#' @return Vector of densities
normaldens <- function(age, mu, lnsigma) {
  sigma <- exp(lnsigma)  # Transform from log scale
  return(dnorm(age, mean = mu, sd = sigma))
}


#' Lognormal Distribution Density
#' 
#' Computes probability density for a lognormal distribution
#' Useful when age distribution is right-skewed
#'
#' @param age Vector of ages (must be positive)
#' @param mu Mean of log(age)
#' @param lnsigma Log of standard deviation of log(age)
#' @return Vector of densities
lognormaldens <- function(age, mu, lnsigma) {
  sigma <- exp(lnsigma)
  # Lognormal: density of log(age) is normal, divide by age for transformation
  return(dnorm(log(age), mean = mu, sd = sigma) / age)
}


#' Weibull Distribution Density
#' 
#' Computes Weibull density (common for mortality modeling)
#'
#' @param age Vector of ages
#' @param lnalpha Log of scale parameter
#' @param lnbeta Log of shape parameter
#' @return Vector of densities
weibulldens <- function(age, lnalpha, lnbeta) {
  alpha <- exp(lnalpha)
  beta <- exp(lnbeta)
  
  # Weibull density formula
  result <- (beta / alpha) * (age / alpha)^(beta - 1) * exp(-(age / alpha)^beta)
  return(result)
}


#' Uniform Distribution Density
#' 
#' Non-informative prior: all ages equally likely in range [a, b]
#'
#' @param age Vector of ages
#' @param a Lower bound
#' @param b Upper bound
#' @return Vector of densities (1/(b-a) within bounds, 0 outside)
uniformdens <- function(age, a, b) {
  result <- ifelse(age >= a & age <= b, 1 / (b - a), 0)
  return(result)
}


#' Gompertz Distribution Density
#' 
#' Classic mortality model where hazard increases exponentially with age
#'
#' @param age Vector of ages
#' @param lnalpha Log of baseline mortality
#' @param lnbeta Log of rate of mortality increase
#' @param aref Reference age for parameterization (default 0)
#' @return Vector of densities
gompertzdens <- function(age, lnalpha, lnbeta, aref = 0) {
  alpha <- exp(lnalpha)
  beta <- exp(lnbeta)
  
  # Reparameterize alpha relative to reference age
  alphatrans <- alpha / exp(beta * aref)
  
  # Compute e^(beta * age)
  ebx <- exp(beta * age)
  
  # Gompertz density
  result <- (alphatrans * ebx) * exp((alphatrans / beta) * (1 - ebx))
  
  # Handle numerical issues (set to 0 if NA)
  result[is.na(result)] <- 0
  
  return(result)
}


#' Gompertz-Makeham Distribution Density
#' 
#' Extension of Gompertz with additional age-independent mortality component
#'
#' @param age Vector of ages
#' @param lnalpha Log of baseline age-dependent mortality
#' @param lnbeta Log of rate of mortality increase
#' @param lnlambda Log of age-independent mortality
#' @param aref Reference age for parameterization (default 0)
#' @return Vector of densities
gompertzmakehamdens <- function(age, lnalpha, lnbeta, lnlambda, aref = 0) {
  alpha <- exp(lnalpha)
  beta <- exp(lnbeta)
  lambda <- exp(lnlambda)
  
  # Reparameterize alpha
  alphatrans <- alpha / exp(beta * aref)
  
  # Compute e^(beta * age)
  ebx <- exp(beta * age)
  
  # Gompertz-Makeham density
  result <- (alphatrans * ebx + lambda) * 
    exp(-(lambda * age) + (alphatrans / beta) * (1 - ebx))
  
  result[is.na(result)] <- 0
  
  return(result)
}


#' Spline Density
#' 
#' Flexible density using cubic splines on log scale
#' Allows for complex, non-parametric age distributions
#'
#' @param age Vector of ages
#' @param knots Vector of knot locations for spline
#' @param beta_params Vector of spline coefficients
#' @return Vector of densities
splinedens <- function(age, knots, beta_params) {
  
  # Create cubic spline basis
  spline_basis <- ns(age, knots = knots, Boundary.knots = range(knots))
  
  # Linear predictor on log scale
  log_density <- as.numeric(spline_basis %*% beta_params)
  
  # Exponentiate to get density
  result <- exp(log_density)
  
  return(result)
}


#' Polygon Density
#' 
#' Piecewise linear density (linear interpolation between specified points)
#' Useful for empirical or semi-parametric priors
#'
#' @param age Vector of ages
#' @param x Vector of x-coordinates (ages) defining polygon vertices
#' @param lnalpha_params Log of y-coordinates (densities) at interior points
#' @param fix Index of fixed point (set to exp(0)=1 for identifiability)
#' @return Vector of densities
polygondens <- function(age, x, lnalpha_params, fix = 2) {
  
  np <- length(x)
  
  # Initialize y-coordinates (densities at each x point)
  y <- numeric(np)
  y[1] <- 0  # Density is 0 at boundaries
  y[np] <- 0
  
  # Fill in interior points
  interior_idx <- 2:(np - 1)
  param_idx <- 1
  for (i in interior_idx) {
    if (i == fix) {
      y[i] <- 1  # Fixed point for identifiability
    } else {
      y[i] <- exp(lnalpha_params[param_idx])
      param_idx <- param_idx + 1
    }
  }
  
  # Linear interpolation between points
  result <- approx(x = x, y = y, xout = age, rule = 2)$y
  
  return(result)
}


# ==============================================================================
# MAIN ESTIMATION FUNCTION
# ==============================================================================

#' Compute f(c|a) - Joint Probability of All Characteristics Given Age
#' 
#' For each individual, computes the probability of observing their complete
#' set of skeletal characteristics given each possible age in the grid.
#' Assumes conditional independence: f(c1,c2,...|a) = f(c1|a) × f(c2|a) × ...
#'
#' @param data Data frame with characteristic variables
#' @param char_vars Character vector of variable names
#' @param char_models List of model specifications for each characteristic
#' @param age Vector of ages to evaluate at
#' @return Matrix: rows = individuals, columns = ages, values = f(c|a)
computefca <- function(data, char_vars, char_models, age) {
  
  n_obs <- nrow(data)
  n_ages <- length(age)
  
  # Initialize result matrix: f(c|a) for each individual at each age
  fca <- matrix(1, nrow = n_obs, ncol = n_ages)
  
  # For each characteristic variable
  for (i in seq_along(char_vars)) {
    var_name <- char_vars[i]
    model_spec <- char_models[[i]]
    
    # Extract model function and parameters
    model_func <- model_spec$func
    model_params <- model_spec$params
    
    # Create expanded data: each individual × each age
    expanded_data <- expand.grid(
      obs_id = 1:n_obs,
      age_id = 1:n_ages
    )
    
    expanded_data$age <- age[expanded_data$age_id]
    expanded_data$char_value <- data[[var_name]][expanded_data$obs_id]
    
    # Compute f(characteristic|age) for all combinations
    expanded_data$prob <- do.call(model_func, c(
      list(vary = expanded_data$char_value, 
           vara = expanded_data$age),
      model_params
    ))
    
    # Reshape back to matrix form
    prob_matrix <- matrix(expanded_data$prob, nrow = n_obs, ncol = n_ages, 
                         byrow = FALSE)
    
    # Multiply into joint probability (conditional independence assumption)
    # Only for non-missing observations
    non_missing <- !is.na(data[[var_name]])
    fca[non_missing, ] <- fca[non_missing, ] * prob_matrix[non_missing, ]
  }
  
  return(fca)
}


#' Population Age Distribution Estimation
#' 
#' Main function: estimates population age distribution f(a) from skeletal data
#' Uses maximum likelihood to find f(a) that best explains observed characteristics
#'
#' The likelihood for individual i is:
#'   L_i = sum over ages a of [f(c_i|a) × f(a)]
#' where f(c_i|a) is computed from characteristic models and f(a) is estimated
#'
#' @param data Data frame with skeletal characteristic observations
#' @param char_vars Character vector of variable names to use
#' @param char_models List of model specifications for each characteristic
#' @param densage_func Function for age density (e.g., normaldens)
#' @param densage_params Named list of starting parameters for age distribution
#' @param gridage Vector of ages to evaluate (e.g., 1:120)
#' @param method Optimization method (default "L-BFGS-B")
#' @param use_trapezoid Use trapezoidal rule for integration (default TRUE)
#' @return List with optimization results and estimated parameters
popest <- function(data, char_vars, char_models, 
                  densage_func, densage_params,
                  gridage, method = "L-BFGS-B",
                  use_trapezoid = TRUE) {
  
  # Remove observations with all characteristics missing
  data_clean <- data %>%
    filter(if_any(all_of(char_vars), ~ !is.na(.)))
  
  # Compute f(c|a) for all individuals at all ages in grid
  message("Computing f(c|a) for all individuals...")
  fca_matrix <- computefca(data_clean, char_vars, char_models, gridage)
  
  # Compute integration weights (trapezoidal rule if requested)
  if (use_trapezoid) {
    n_ages <- length(gridage)
    width <- numeric(n_ages)
    
    # Interior points: average of spacing on left and right
    for (i in 2:(n_ages - 1)) {
      width[i] <- (gridage[i + 1] - gridage[i - 1]) / 2
    }
    
    # Boundary points: half spacing to neighbor
    width[1] <- (gridage[2] - gridage[1]) / 2
    width[n_ages] <- (gridage[n_ages] - gridage[n_ages - 1]) / 2
  } else {
    # Uniform weights (simple sum)
    width <- rep(1, length(gridage))
  }
  
  # Define negative log-likelihood function to minimize
  neg_loglik <- function(params) {
    
    # Convert parameter vector to named list
    param_list <- as.list(params)
    names(param_list) <- names(densage_params)
    
    # Compute f(a) for all ages in grid using current parameters
    fa <- do.call(densage_func, c(list(age = gridage), param_list))
    
    # Apply integration weights
    fa_weighted <- fa * width
    
    # Normalize so f(a) integrates to 1
    fa_normalized <- fa_weighted / sum(fa_weighted)
    
    # For each individual, compute likelihood:
    # L_i = sum_a [f(c_i|a) × f(a)]
    likelihood <- as.numeric(fca_matrix %*% fa_normalized) # Matrix multiplication to calculate the marginal likelihood for each individual
    ## "What's the probability of seeing individual i's skeletal characteristics, 
    ## averaging over all possible ages they could have been, weighted by how common each age is in the population?"
    
    # Log-likelihood
    loglik <- sum(log(likelihood))
    
    # Return negative (for minimization)
    return(-loglik)
  }
  
  # Starting values
  start_params <- unlist(densage_params)
  
  # Optimize
  message("Optimizing age distribution parameters...")
  opt_result <- optim(
    par = start_params,
    fn = neg_loglik,
    method = method,
    control = list(trace = 1, maxit = 1000)
  )
  
  # Extract results
  estimated_params <- as.list(opt_result$par)
  names(estimated_params) <- names(densage_params)
  
  result <- list(
    params = estimated_params,
    loglik = -opt_result$value,
    convergence = opt_result$convergence,
    densage_func = densage_func,
    gridage = gridage,
    data = data_clean,
    char_vars = char_vars,
    char_models = char_models,
    fca_matrix = fca_matrix
  )
  
  class(result) <- "popest"
  
  return(result)
}


# ==============================================================================
# POST-ESTIMATION FUNCTIONS
# ==============================================================================

#' Describe Age Distribution
#' 
#' Computes summary statistics of estimated age distribution f(a)
#' Such as mean, median, percentiles, probabilities of age ranges
#'
#' @param popest_result Result object from popest()
#' @param stats Character vector of statistics to compute
#'   Options: "mean", "median", "sd", "p10", "p90", etc.
#'   Or probabilities like "probgt(40)" for P(age > 40)
#' @param plot Logical, whether to create visualization
#' @return Data frame of computed statistics (and plot if requested)
describedist <- function(popest_result, stats = c("mean", "median", "sd"), 
                        plot = FALSE) {
  
  # Extract components
  gridage <- popest_result$gridage
  densage_func <- popest_result$densage_func
  params <- popest_result$params
  
  # Compute f(a) at grid points
  fa <- do.call(densage_func, c(list(age = gridage), params))
  
  # Normalize
  fa <- fa / sum(fa)
  
  # Create data frame for computation
  age_dist <- data.frame(age = gridage, prob = fa)
  
  # Compute requested statistics
  results <- list()
  
  for (stat in stats) {
    if (stat == "mean") {
      results$mean <- sum(age_dist$age * age_dist$prob)
      
    } else if (stat == "median") {
      cum_prob <- cumsum(age_dist$prob)
      results$median <- age_dist$age[which(cum_prob >= 0.5)[1]]
      
    } else if (stat == "sd") {
      mean_age <- sum(age_dist$age * age_dist$prob)
      results$sd <- sqrt(sum((age_dist$age - mean_age)^2 * age_dist$prob))
      
    } else if (grepl("^p[0-9]+$", stat)) {
      # Percentile (e.g., "p10" for 10th percentile)
      pct <- as.numeric(gsub("p", "", stat)) / 100
      cum_prob <- cumsum(age_dist$prob)
      results[[stat]] <- age_dist$age[which(cum_prob >= pct)[1]]
      
    } else if (grepl("^probgt\\([0-9]+\\)$", stat)) {
      # Probability greater than threshold
      threshold <- as.numeric(gsub("probgt\\(|\\)", "", stat))
      results[[stat]] <- sum(age_dist$prob[age_dist$age > threshold]) * 100
      
    } else if (grepl("^probge\\([0-9]+\\)$", stat)) {
      # Probability greater than or equal to threshold
      threshold <- as.numeric(gsub("probge\\(|\\)", "", stat))
      results[[stat]] <- sum(age_dist$prob[age_dist$age >= threshold]) * 100
      
    } else if (grepl("^problt\\([0-9]+\\)$", stat)) {
      # Probability less than threshold
      threshold <- as.numeric(gsub("problt\\(|\\)", "", stat))
      results[[stat]] <- sum(age_dist$prob[age_dist$age < threshold]) * 100
      
    } else if (grepl("^proble\\([0-9]+\\)$", stat)) {
      # Probability less than or equal to threshold
      threshold <- as.numeric(gsub("proble\\(|\\)", "", stat))
      results[[stat]] <- sum(age_dist$prob[age_dist$age <= threshold]) * 100
    }
  }
  
  # Convert to data frame
  results_df <- as.data.frame(results)
  
  # Create plot if requested
  if (plot) {
    p <- ggplot(age_dist, aes(x = age, y = prob * 100)) +
      geom_line() +
      labs(x = "Age", y = "% per year",
           title = "Estimated Population Age Distribution") +
      theme_minimal()
    
    print(p)
  }
  
  return(results_df)
}


#' Individual Age Estimates
#' 
#' Computes posterior age distribution for each individual given their 
#' skeletal characteristics: f(a|c_i) ∝ f(c_i|a) × f(a)
#'
#' Returns posterior mean, SD, and/or mode for each individual
#'
#' @param popest_result Result object from popest()
#' @param compute_mean Compute posterior mean age?
#' @param compute_sd Compute posterior standard deviation?
#' @param compute_mode Compute posterior mode (most likely age)?
#' @return Data frame with individual-level estimates
indest <- function(popest_result, 
                  compute_mean = TRUE,
                  compute_sd = FALSE, 
                  compute_mode = FALSE) {
  
  # Extract components
  gridage <- popest_result$gridage
  densage_func <- popest_result$densage_func
  params <- popest_result$params
  fca_matrix <- popest_result$fca_matrix
  n_obs <- nrow(fca_matrix)
  
  # Compute prior f(a)
  fa <- do.call(densage_func, c(list(age = gridage), params))
  fa <- fa / sum(fa)  # Normalize
  
  # Compute posterior for each individual: f(a|c_i) ∝ f(c_i|a) × f(a)
  posterior <- fca_matrix * matrix(fa, nrow = n_obs, ncol = length(gridage), 
                                   byrow = TRUE)
  
  # Normalize each posterior distribution
  posterior_sums <- rowSums(posterior)
  posterior <- posterior / posterior_sums
  
  # Initialize result data frame
  results <- data.frame(individual = 1:n_obs)
  
  # Compute posterior mean
  if (compute_mean) {
    # E[a|c] = sum_a [a × f(a|c)]
    results$postmean <- as.numeric(posterior %*% gridage)
  }
  
  # Compute posterior standard deviation
  if (compute_sd) {
    if (!compute_mean) {
      postmean <- as.numeric(posterior %*% gridage)
    } else {
      postmean <- results$postmean
    }
    
    # SD[a|c] = sqrt(E[(a - E[a|c])^2 | c])
    age_matrix <- matrix(gridage, nrow = n_obs, ncol = length(gridage), 
                        byrow = TRUE)
    squared_dev <- (age_matrix - postmean)^2
    results$postsd <- sqrt(rowSums(posterior * squared_dev))
  }
  
  # Compute posterior mode
  if (compute_mode) {
    # Most likely age = age with highest posterior probability
    max_idx <- apply(posterior, 1, which.max)
    results$postmode <- gridage[max_idx]
  }
  
  return(results)
}


# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================

#' Example: Simulated Skeletal Data
#' 
#' This example demonstrates the workflow:
#' 1. Define models for skeletal characteristics f(y|a)
#' 2. Estimate population age distribution f(a)
#' 3. Compute individual age estimates
example_usage <- function() {
  
  # Simulate example data
  set.seed(123)
  n <- 100
  
  # True ages (unknown in practice)
  true_ages <- sample(20:80, n, replace = TRUE)
  
  # Simulate three ordinal skeletal characteristics
  # Each follows ordinal regression model with some error
  
  # Characteristic 1: stages 1-7
  y1 <- pmax(1, pmin(7, round((true_ages - 20) / 10) + sample(1:7, n, replace = TRUE, 
                                                               prob = c(0.3, 0.3, 0.2, 0.1, 0.05, 0.03, 0.02))))
  
  # Characteristic 2: stages 1-7  
  y2 <- pmax(1, pmin(7, round((true_ages - 15) / 10) + sample(1:7, n, replace = TRUE,
                                                               prob = c(0.25, 0.25, 0.25, 0.15, 0.05, 0.03, 0.02))))
  
  # Characteristic 3: stages 1-7
  y3 <- pmax(1, pmin(7, round((true_ages - 25) / 10) + sample(1:7, n, replace = TRUE,
                                                               prob = c(0.3, 0.25, 0.2, 0.15, 0.05, 0.03, 0.02))))
  
  # Create data frame
  data <- data.frame(y1, y2, y3)
  
  # Define models for each characteristic
  # Using ordinal logistic regression model
  char_models <- list(
    list(func = oregressprob,
         params = list(distf = "logistic", from = 1, to = 7,
                      kappa = c(25, 35, 45, 55, 65, 75),
                      sigma = 10.0)),
    
    list(func = oregressprob,
         params = list(distf = "logistic", from = 1, to = 7,
                      kappa = c(25, 35, 45, 55, 65, 75),
                      sigma = 10.0)),
    
    list(func = oregressprob,
         params = list(distf = "logistic", from = 1, to = 7,
                      kappa = c(25, 35, 45, 55, 65, 75),
                      sigma = 10.0))
  )
  
  # Estimate population age distribution
  # Using lognormal prior with initial guesses
  result <- popest(
    data = data,
    char_vars = c("y1", "y2", "y3"),
    char_models = char_models,
    densage_func = lognormaldens,
    densage_params = list(mu = 3.5, lnsigma = -1),
    gridage = seq(1, 120, by = 1)
  )
  
  # Describe estimated age distribution
  summary_stats <- describedist(
    result,
    stats = c("mean", "median", "sd", "p10", "p90", 
             "probgt(40)", "problt(50)"),
    plot = TRUE
  )
  
  print("Population age distribution summary:")
  print(summary_stats)
  
  # Compute individual age estimates
  individual_estimates <- indest(
    result,
    compute_mean = TRUE,
    compute_sd = TRUE,
    compute_mode = TRUE
  )
  
  print("First 10 individual age estimates:")
  print(head(individual_estimates, 10))
  
  # Compare to true ages (only possible in simulation)
  comparison <- data.frame(
    true_age = true_ages,
    estimated_age = individual_estimates$postmean,
    post_sd = individual_estimates$postsd
  )
  
  print("Estimation accuracy:")
  print(paste("RMSE:", round(sqrt(mean((comparison$true_age - comparison$estimated_age)^2)), 2)))
  print(paste("MAE:", round(mean(abs(comparison$true_age - comparison$estimated_age)), 2)))
  
  return(list(result = result, individual_estimates = individual_estimates,
              comparison = comparison))
}

# ==============================================================================
# END OF TOOLKIT
# ==============================================================================
