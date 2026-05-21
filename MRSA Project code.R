# Load libraries
library(deSolve)
library(readxl)
library(tidyverse)
library(ggplot2)

#---------------------------
# SCIR Model Definition
#---------------------------
SCIR <- function(time, current_state, params){
  with(as.list(c(current_state, params)), {
    N <- S + C + I + R
    dS <- p * N - (lambda_C * C * S) / N - (lambda_I * C * S) / N - 
      (beta_C * I * S) / N - (beta_I * I * S) / N + alpha * C + 
      mu * R - delta * S
    dC <- (lambda_C * C * S) / N + (lambda_I * C * S) / N - 
      alpha * C - gamma * C - delta * C
    dI <- gamma * C + (beta_C * I * S) / N + (beta_I * I * S) / N - 
      nu * I - delta * I
    dR <- nu * I - mu * R - delta * R
    return(list(c(dS, dC, dI, dR)))
  })
}

#---------------------------
# Load & Prepare Data
#---------------------------
infections <- read_xls("~/Dataset_S1.xls")
infections$month <- as.Date(paste0(infections$begin_year, "-", infections$begin_month, "-01"))

monthly_cases_train <- infections %>%
  filter(begin_year %in% 2004:2006) %>%
  group_by(month) %>%
  summarise(I_obs = n()) %>%
  arrange(month)

monthly_cases_train$time_index <- as.integer(difftime(monthly_cases_train$month, 
                                                      min(monthly_cases_train$month), 
                                                      units = "days") / 30)

times <- monthly_cases_train$time_index
time_days <- seq(0, max(monthly_cases_train$time_index * 30), by = 1)
I_obs <- monthly_cases_train$I_obs

#---------------------------
# Model Simulator
#---------------------------
run_model <- function(params, times, init_state) {
  out <- ode(y = init_state, times = time_days, func = SCIR, parms = params)
  as.data.frame(out)
}

#---------------------------
# Log-Likelihood (Poisson)
#---------------------------
log_likelihood <- function(params, time_days, init_state, I_obs, time_index_obs) {
  model_out <- run_model(params, time_days, init_state)
  model_out$time_index <- floor(model_out$time / 30)
  
  monthly_pred <- model_out %>%
    group_by(time_index) %>%
    summarise(I_pred = sum(I)) %>%
    filter(time_index %in% time_index_obs)
  
  if (nrow(monthly_pred) != length(I_obs)) return(-Inf)
  
  # Avoid invalid values
  eps <- 1e-10
  monthly_pred$I_pred <- pmax(monthly_pred$I_pred, eps)
  
  ll <- sum(dpois(I_obs, lambda = monthly_pred$I_pred, log = TRUE))
  if (is.nan(ll)) return(-Inf)
  
  return(ll)
}

#---------------------------
# MCMC Sampler
#---------------------------
mcmc <- function(init_params, init_state, times, I_obs, time_index_obs, n_iter = 10000, proposal_sd = 0.5) {
  chain <- matrix(NA, nrow = n_iter, ncol = length(init_params))
  colnames(chain) <- names(init_params)
  current_params <- init_params
  current_ll <- log_likelihood(current_params, times, init_state, I_obs, time_index_obs)
  
  for (i in 1:n_iter) {
    proposed_params <- abs(current_params + rnorm(length(init_params), mean = 0, sd = proposal_sd))
    names(proposed_params) <- names(init_params)
    
    proposed_ll <- log_likelihood(proposed_params, times, init_state, I_obs, time_index_obs)
    log_accept_ratio <- proposed_ll - current_ll
    
    if (!is.nan(log_accept_ratio) && log(runif(1)) < log_accept_ratio) {
      current_params <- proposed_params
      current_ll <- proposed_ll
    }
    
    chain[i, ] <- current_params
  }
  
  return(as.data.frame(chain))
}

#---------------------------
# Initial Parameters
#---------------------------
init_state <- c(S = 999, C = 1, I = 0, R = 0)
init_params <- c(lambda_C = 1.46e-2, lambda_I = 1.76e-2, 
                 beta_C = 7.30e-6, beta_I = 1.49e-2, 
                 alpha = 9.80e-3, gamma = 1.20e-4, 
                 nu = 14.60, mu = 7.70, delta = 0.01, p = 0.01)

#---------------------------
# Run MCMC
#---------------------------
set.seed(112)
chain <- mcmc(init_params, init_state, time_days, I_obs, time_index_obs = monthly_cases_train$time_index)

# Get posterior means
posterior_means <- colMeans(chain)

#---------------------------
# Simulate with Posterior
#---------------------------
sim_result <- run_model(posterior_means, time_days, init_state)
sim_result$time_index <- floor(sim_result$time / 30)

monthly_pred <- sim_result %>%
  group_by(time_index) %>%
  summarise(I_pred = sum(I)) %>%
  filter(time_index %in% monthly_cases_train$time_index)

comparison <- merge(monthly_cases_train, monthly_pred, by = "time_index")

#---------------------------
# Plot Results
#---------------------------
plot(comparison$time_index, comparison$I_obs, type = 'p', col = 'red', pch = 16,
     ylab = "Number of Infections", xlab = "Month Index", main = "Observed vs Predicted Infections (2004-2006)")
lines(comparison$time_index, comparison$I_pred, col = 'blue', lwd = 2)
legend("topleft", legend = c("Observed", "Predicted"), col = c("red", "blue"), lwd = c(NA,2), pch = c(16, NA))

#---------------------------
# Simulate Next 24 Months
#---------------------------
#---------------------------
# Prepare Test Data
#---------------------------
monthly_cases_test <- infections %>%
  filter(begin_year %in% 2007:2008) %>%
  group_by(month) %>%
  summarise(I_obs = n()) %>%
  arrange(month)

# Align test time index with training
monthly_cases_test$time_index <- as.integer(difftime(monthly_cases_test$month,
                                                     min(monthly_cases_train$month),
                                                     units = "days") / 30)

#---------------------------
# Final State from Training Simulation
#---------------------------
final_training_state <- sim_result %>%
  filter(time == max(time)) %>%
  select(S, C, I, R) %>%
  unlist()

# Use this as initial state for test sim
init_state_test <- c(S = 1400, C = 851, I = 45, R = 86)

#---------------------------
# Define Test Simulation Time
#---------------------------
test_duration_months <- 24
test_duration_days <- test_duration_months * 30
time_days_test <- seq(0, test_duration_days, by = 1)
expected_time_indices <- 0:(test_duration_days %/% 30)
n_timepoints <- length(expected_time_indices)

#---------------------------
# Sample 500 Posterior Parameter Sets
#---------------------------
n_samples <- 500
posterior_samples <- chain[sample(1:nrow(chain), n_samples, replace = TRUE), ]

#---------------------------
# Run Simulations
#---------------------------
I_predictions_test <- matrix(NA, nrow = n_samples, ncol = n_timepoints)

for (i in 1:n_samples) {
  params <- as.numeric(posterior_samples[i, ])
  names(params) <- names(posterior_means)
  
  sim <- run_model(params, time_days_test, init_state_test)
  sim$time_index <- floor(sim$time / 30)
  
  monthly_I <- sim %>%
    group_by(time_index) %>%
    summarise(I_pred = sum(I)) %>%
    right_join(data.frame(time_index = expected_time_indices), by = "time_index") %>%
    arrange(time_index) %>%
    mutate(I_pred = replace_na(I_pred, 0))
  
  I_predictions_test[i, ] <- monthly_I$I_pred
}

#---------------------------
# Filter Simulations: Max Infections ≤ 2000
#---------------------------
valid_rows <- apply(I_predictions_test, 1, max) >= 0 & apply(I_predictions_test, 1, max) <= 2000
I_predictions_filtered <- I_predictions_test[valid_rows, ]

# Report how many were kept
cat("Simulations kept:", sum(valid_rows), "/", n_samples, "\n")

#---------------------------
# Posterior Predictive Summaries
#---------------------------
summary_test_pred_filtered <- data.frame(
  time_index = expected_time_indices,
  mean = apply(I_predictions_filtered, 2, mean),
  lower_95 = apply(I_predictions_filtered, 2, quantile, probs = 0.025),
  upper_95 = apply(I_predictions_filtered, 2, quantile, probs = 0.975)
)

#---------------------------
# Prepare Test Observations
#---------------------------
monthly_cases_test <- infections %>%
  filter(begin_year %in% 2007:2008) %>%
  group_by(month) %>%
  summarise(I_obs = n()) %>%
  arrange(month)

monthly_cases_test$time_index <- as.integer(difftime(monthly_cases_test$month,
                                                     min(monthly_cases_train$month),
                                                     units = "days") / 30)

# Align test indices to start at 0
monthly_cases_test_shifted <- monthly_cases_test %>%
  mutate(time_index = time_index - min(time_index))

#---------------------------
# Merge and Plot
#---------------------------
comparison_test <- merge(monthly_cases_test_shifted, summary_test_pred_filtered, by = "time_index")

plot(comparison_test$time_index, comparison_test$I_obs, type = "p", col = "red", pch = 16,
     ylab = "Number of Infections", xlab = "Month Index",
     main = "Observed vs Predicted Infections (2007-2008)",
     ylim = c(0, 2000))
lines(comparison_test$time_index, comparison_test$mean, col = "blue", lwd = 2)
lines(comparison_test$time_index, comparison_test$lower_95, col = "gray60", lty = 2)
lines(comparison_test$time_index, comparison_test$upper_95, col = "gray60", lty = 2)
polygon(c(comparison_test$time_index, rev(comparison_test$time_index)),
        c(comparison_test$upper_95, rev(comparison_test$lower_95)),
        col = rgb(0.6, 0.6, 0.6, 0.3), border = NA)
legend("bottomright", legend = c("Observed", "Mean Prediction", "95% CI"),
       col = c("red", "blue", "gray60"), pch = c(16, NA, NA), lty = c(NA, 1, 2), lwd = c(NA, 2, 1))

#---------------------------
# Evaluation Metrics
#---------------------------
RMSE_test <- sqrt(mean((comparison_test$I_obs - comparison_test$mean)^2))
MAE_test <- mean(abs(comparison_test$I_obs - comparison_test$mean))
cat("Test RMSE:", RMSE_test, "\nTest MAE:", MAE_test, "\n")

library(ggplot2)
library(tidyr)

# Add an iteration column to the chain data
chain$iteration <- 1:nrow(chain)

# Reshape the MCMC chain to long format
chain_long <- pivot_longer(chain, 
                           cols = -iteration,  # Keep 'iteration' as it is
                           names_to = "parameter", 
                           values_to = "value")

# Plot the trace plots for each parameter
ggplot(chain_long, aes(x = iteration, y = value)) + 
  geom_line() + 
  facet_wrap(~parameter, scales = "free_y") + 
  theme_minimal() +
  labs(title = "MCMC Trace Plots", x = "Iteration", y = "Parameter Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
library(ggplot2)
library(tidyr)

# Add an iteration column to the chain data
chain$iteration <- 1:nrow(chain)

# Reshape the MCMC chain to long format
chain_long <- pivot_longer(chain, 
                           cols = -iteration,  # Keep 'iteration' as it is
                           names_to = "parameter", 
                           values_to = "value")

# Plot the trace plots for each parameter
ggplot(chain_long, aes(x = iteration, y = value)) + 
  geom_line() + 
  facet_wrap(~parameter, scales = "free_y") + 
  theme_minimal() +
  labs(title = "MCMC Trace Plots", x = "Iteration", y = "Parameter Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
