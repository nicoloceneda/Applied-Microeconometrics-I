# ------------------------------------------------------------------------------
# Prepare the environment
# ------------------------------------------------------------------------------

# Load the libraries

library(causalForest)
library(mgcv)
library(randomForestCI)
library(FNN)
library(Hmisc)
library(xtable)

# Clear the memory

rm(list = ls())

# ------------------------------------------------------------------------------
# Variables 
# ------------------------------------------------------------------------------

# Define the variables

n = 500                                                                         # Number of samples {n}
n.test = 100                                                                    # Number of test samples
dvals = c(2, 5, 10, 15, 20, 30)                                                 # Number of features {d}
sigma = 1                                                                       # Volatility {sigma}

ntree = 100 #1000                                                               # Number of trees {B}
simu.reps = 100 #500                                                            # Number of simulations

# ------------------------------------------------------------------------------
# Simulation function:
# ------------------------------------------------------------------------------

# Define the functions in (27)

propensity = function(x) {                                                      # Function for treatment propensity {e(x)}
  0.25 * (1 + dbeta(x[1], 2, 4))
}

baseline = function(x) {                                                        # Function for main effect {m(x)}
  2 * x[1] - 1
}

simu.fun = function(d) {
  
  # ----------------------------------------------------------------------------
  # Generate data
  # ----------------------------------------------------------------------------
  
  # Generate X
  X = matrix(runif(n * d, 0, 1), n, d)                                          # Matrix of features {X ~ U([0,1])}
  X.test = matrix(runif(n.test * d, 0, 1), n.test, d)                           # Matrix of test features {X_test ~ U([0,1])}
  
  # Generate Y
  Y = apply(X, 1, baseline) + sigma * rnorm(n)                                  # Vector of outcomes {Y ~ N(E[Y|X],1)}
  
  # Generate W
  e = apply(X, 1, propensity)                                                   # Treatment propensity in (27)
  W = rbinom(n, 1, e)                                                           # Treatment assignment with Bernoulli distribution
  
  # ----------------------------------------------------------------------------
  # Random forest
  # ----------------------------------------------------------------------------
  
  # Create a random forest with propensity trees [F1]
  forest = propensityForest(X, Y, W, num.trees=ntree, sample.size=n/10, nodesize=1)                                      
  
  # Compute predictions on the test set
  predictions = predict(forest, X.test)
  
  # Compute the confidence interval on the test set
  forest.ci = randomForestInfJack(forest, X.test, calibrate = TRUE)
  
  # Evaluate performance of random forests
  se.hat = sqrt(forest.ci$var.hat)                                              # Standard error
  rf.covered = mean(abs(predictions) <= 1.96 * se.hat)                          # Expected coverage (mean of predictions contained in CI)
  rf.mse = mean(predictions^2)                                                  # Expected mean-squared error (true treatment effect tau(x) = 0)
  
  # ----------------------------------------------------------------------------
  # K nearest neighbours (10 neigh)
  # ----------------------------------------------------------------------------
  
  # Compute predictions on the test
  k.small = 10                                                                  # Number of neighbours considered
  knn.0.mu = knn.reg(X[W==0,], X.test, Y[W==0], k = k.small)$pred               # Predicted values for E[Y|W=0]       
  knn.1.mu = knn.reg(X[W==1,], X.test, Y[W==1], k = k.small)$pred               # Predicted values for E[Y|W=1]
  knn.tau = knn.1.mu - knn.0.mu                                                 # tau(x) = E[Y|W=1] - E[Y|W=0]
  
  knn.0.mu2 = knn.reg(X[W==0,], X.test, Y[W==0]^2, k = k.small)$pred            # Predicted values for E[Y^2|W=0] 
  knn.1.mu2 = knn.reg(X[W==1,], X.test, Y[W==1]^2, k = k.small)$pred            # Predicted values for E[Y^2|W=1]
  
  # Evaluate performance of k nearest neighbor
  knn.0.var = (knn.0.mu2 - knn.0.mu^2) / (k.small - 1)                          # Var[Y|W=0]
  knn.1.var = (knn.1.mu2 - knn.1.mu^2) / (k.small - 1)                          # Var[Y|W=0]
  knn.se = sqrt(knn.0.var + knn.1.var)                                          # Standard error
  
  knn.covered.v1 = mean(abs(knn.tau) <= 1.96 * knn.se)                          # Expected coverage (mean of predictions contained in CI)
  knn.covered.v2 = mean(abs(knn.tau) <= 1.96 * sqrt(mean(knn.0.var+knn.1.var))) # Expected coverage (mean of predictions contained in CI)
  
  knn.mse = mean(knn.tau^2)                                                     # Expected mean-squared error (true treatment effect tau(x) = 0)
  
  # ----------------------------------------------------------------------------
  # K nearest neighbours (100 neigh)
  # ----------------------------------------------------------------------------
  
  # Compute predictions on the test
  k.big = 100                                                                   # Number of neighbours considered
  knnbig.0.mu = knn.reg(X[W==0,], X.test, Y[W==0], k = k.big)$pred              # Predicted values for E[Y|W=0] 
  knnbig.1.mu = knn.reg(X[W==1,], X.test, Y[W==1], k = k.big)$pred              # Predicted values for E[Y|W=1]
  knnbig.tau = knnbig.1.mu - knnbig.0.mu                                        # tau(x) = E[Y|W=1] - E[Y|W=0]
  
  knnbig.0.mu2 = knn.reg(X[W==0,], X.test, Y[W==0]^2, k = k.big)$pred           # Predicted values for E[Y^2|W=0]
  knnbig.1.mu2 = knn.reg(X[W==1,], X.test, Y[W==1]^2, k = k.big)$pred           # Predicted values for E[Y^2|W=1]
  
  # Evaluate performance of k nearest neighbor
  knnbig.0.var = (knnbig.0.mu2 - knnbig.0.mu^2) / (k.big - 1)                   # Var[Y|W=0]
  knnbig.1.var = (knnbig.1.mu2 - knnbig.1.mu^2) / (k.big - 1)                   # Var[Y|W=0]
  knnbig.se = sqrt(knnbig.0.var + knnbig.1.var)                                 # Standard error
  
  knnbig.covered.v1 = mean(abs(knnbig.tau) <= 1.96 * knnbig.se)                 # Expected coverage (mean of predictions contained in CI)
  knnbig.covered.v2=mean(abs(knnbig.tau)<=1.96*sqrt(mean(knnbig.0.var+knnbig.1.var))) # Expected coverage (mean of predictions contained in CI)
                                                               
  knnbig.mse = mean(knnbig.tau^2)                                               # Expected mean-squared error (true treatment effect tau(x) = 0)
  
  # ----------------------------------------------------------------------------
  # Results
  # ----------------------------------------------------------------------------
  
  # Save the results
  c(rf.covered = rf.covered,                                                    # RT: expected coverage                     
    rf.mse = rf.mse,                                                            # RT: expected mean-squared error 
    knn.covered.v1 = knn.covered.v1,                                            # KNN_10: expected coverage
    knn.covered.v2 = knn.covered.v2,                                            # KNN_10: expected coverage
    knn.mse = knn.mse,                                                          # KNN_10: expected mean-squared error
    knnbig.covered.v1 = knnbig.covered.v1,                                      # KNN_100: expected coverage
    knnbig.covered.v2 = knnbig.covered.v2,                                      # KNN_100: expected coverage
    knnbig.mse = knnbig.mse)                                                    # KNN_100: expected mean-squared error
}

# ------------------------------------------------------------------------------
# Simulation run:
# ------------------------------------------------------------------------------

# Run the simulations simulations

results.raw = lapply(dvals, function(d) {
  
  print(paste("NOW RUNNING:", d))                                             
  res.d = sapply(1:simu.reps, function(iter) simu.fun(d))                       # Run {simu.reps} simulations for each {dvals}
  res.fixed = data.frame(t(res.d))                                              # Save the result of the simuls as a dataframe
  
  print(paste("RESULT AT", d, "IS", colMeans(res.fixed)))                       # Print the mean of results
  res.fixed                                                                     # Add the dataframe of each d to results.raw
})

# Compute aggregate results of all simulations

results.condensed = lapply(results.raw, function(RR) {
  RR.mu = colMeans(RR)                                                          # Mean of results
  RR.var = sapply(RR, var) / (nrow(RR) - 1)                                     # Variance of results
  rbind("mu"=RR.mu, "se"=sqrt(RR.var))                                          # Group the results
})

# Display the results

results.condensed                                                               # Print the results
save.image("table1_propensity.RData")                                           # Save the results

# Table of results

results.parsed = lapply(results.condensed, function(RR) {
  apply(RR, 2, function(arg) {
    paste0(round(arg[1], 2), " (", round(100 * arg[2], 0), ")")
  })
})

results.table = data.frame(cbind(d=dvals, Reduce(rbind, results.parsed)))       # Put the result into a dataframe
results.table = results.table[,c(1, 3, 6, 9, 2, 4, 7)]                          # Select the columns
xtab = xtable(results.table)                                                    # Convert R object to xtable object
print(xtab, include.rownames = FALSE)                                           # Print table