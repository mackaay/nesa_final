# Load the dataset
# Ensure 'aluminium.csv' is in your working directory
df <- read.csv("./nesa_test/aluminium.csv")
plot(df$stress ~ df$strain)
x <- df$strain
y <- df$stress
n <- length(y)

# Prior Hyperparameters
beta_0_val <- 30000
v_0_val <- 20000^2
beta0 <- c(beta_0_val, beta_0_val)
V0 <- diag(c(v_0_val, v_0_val))
invV0 <- solve(V0)
a0 <- 4
b0 <- 13.12

# Candidate breakpoints: unique internal strain values
# We exclude the first and last few points to ensure segments are identifiable
taus <- sort(unique(x))
taus <- taus[taus > min(x) & taus < max(x)]

# Vector to store log marginal likelihood for each tau
log_marg_lik <- numeric(length(taus))

for (i in 1:length(taus)) {
  tau <- taus[i]
  
  # Define the design matrix Z based on the piecewise structure:
  # mu = beta1 * min(x, tau) + beta2 * max(0, x - tau)
  z1 <- pmin(x, tau)
  z2 <- pmax(0, x - tau)
  Z <- cbind(z1, z2)
  
  # Update parameters for the Normal-Inverse-Gamma posterior
  Vn <- solve(invV0 + t(Z) %*% Z)
  betan <- Vn %*% (invV0 %*% beta0 + t(Z) %*% y)
  an <- a0 + n/2
  
  # b_n calculation for the conjugate NIG prior
  # bn = b0 + 0.5 * (y'y + beta0' V0^-1 beta0 - betan' Vn^-1 betan)
  term_y <- t(y) %*% y
  term_prior <- t(beta0) %*% invV0 %*% beta0
  term_post <- t(betan) %*% solve(Vn) %*% betan
  bn <- as.numeric(b0 + 0.5 * (term_y + term_prior - term_post))
  
  # Calculate the log marginal likelihood log p(y | tau)
  # Dropping constants that do not depend on tau
  log_marg_lik[i] <- -an * log(bn) + 0.5 * determinant(Vn, logarithm = TRUE)$modulus
}

# Calculate the posterior probabilities of tau (discrete)
# Subtracting max(log_marg_lik) prevents numerical underflow (the 'log-sum-exp' trick)
post_tau_unnorm <- exp(log_marg_lik - max(log_marg_lik))
post_tau <- post_tau_unnorm / sum(post_tau_unnorm)

# Point Estimates
yield_point_mean <- sum(taus * post_tau)
yield_point_map <- taus[which.max(post_tau)]

# Display Results
cat("--- Bayesian Piecewise Linear Model Results ---\n")
cat("Posterior Mean (Expected Yield Point):", yield_point_mean, "\n")
cat("MAP Estimate (Mode of Posterior):    ", yield_point_map, "\n")

# Visualization
plot(taus, post_tau, type = "h", col = "blue", lwd = 2,
     main = "Posterior Distribution of Yield Point (tau)",
     xlab = "Strain (tau)", ylab = "Posterior Probability")
abline(v = yield_point_mean, col = "red", lty = 2, lwd = 2)
legend("topright", legend = c("Posterior Prob", "Posterior Mean"), 
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)


plot(taus, post_tau, type = "h", col = "blue", lwd = 2,
     main = "Posterior Distribution of Yield Point (tau)",
     xlab = "Strain (tau)", ylab = "Posterior Probability",
     xlim = c(0.003, 0.005) )






# Storage for results
log_marg_lik <- numeric(length(taus))
beta_means <- matrix(0, nrow = length(taus), ncol = 2)

# 4. Bayesian Update Loop
for (i in seq_along(taus)) {
  tau <- taus[i]
  
  # Basis functions for the continuous piecewise model
  z1 <- pmin(x, tau)
  z2 <- pmax(0, x - tau)
  Z <- cbind(z1, z2)
  
  # Posterior parameters for beta and sigma^2 given tau
  Vn <- solve(invV0 + t(Z) %*% Z)
  betan <- Vn %*% (invV0 %*% beta0 + t(Z) %*% y)
  an <- a0 + n/2
  
  # bn calculation (sum of squares term for NIG prior)
  term_ss <- as.numeric(t(y) %*% y + t(beta0) %*% invV0 %*% beta0 - t(betan) %*% solve(Vn) %*% betan)
  bn <- b0 + 0.5 * term_ss
  
  # Log-marginal likelihood: log p(y | tau)
  # Includes determinant and bn terms (constants independent of tau can be ignored)
  log_marg_lik[i] <- 0.5 * determinant(Vn, logarithm = TRUE)$modulus - an * log(bn)
  
  # Store the posterior mean of beta for this tau
  beta_means[i, ] <- betan
}
# 5. Compute Posterior Probabilities of tau
post_probs <- exp(log_marg_lik - max(log_marg_lik))
post_probs <- post_probs / sum(post_probs)

# 6. Calculate the Posterior Mean Curve
# We evaluate the curve on a fine grid for smooth plotting
x_grid <- seq(min(x), max(x), length.out = 200)
y_post_mean <- numeric(length(x_grid))

for (i in seq_along(taus)) {
  tau_i <- taus[i]
  prob_i <- post_probs[i]
  b_i <- beta_means[i, ]
  
  # Contribution of this tau to the expected curve:
  # mu(x|tau) = beta1 * min(x, tau) + beta2 * max(0, x - tau)
  y_post_mean <- y_post_mean + prob_i * (b_i[1] * pmin(x_grid, tau_i) + b_i[2] * pmax(0, x_grid - tau_i))
}

# 7. Visualization
plot(x, y, pch = 16, col = rgb(0.2, 0.2, 0.2, 0.5), 
     main = "Aluminium Stress-Strain: Bayesian Piecewise Fit",
     xlab = "Strain", ylab = "Stress")

# Add the posterior mean curve
lines(x_grid, y_post_mean, col = "red", lwd = 3)

# Optional: Add a vertical line for the MAP estimate of the yield point
map_tau <- taus[which.max(post_probs)]
abline(v = map_tau, col = "blue", lty = 2)

legend("bottomright", legend = c("Observed Data", "Posterior Mean Curve", "MAP Yield Point"),
       col = c("gray", "red", "blue"), pch = c(16, NA, NA), lty = c(NA, 1, 2), lwd = c(NA, 3, 1))
