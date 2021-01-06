# ------------------------------------------------------------------------------
# Propensity Forest initialization
# ------------------------------------------------------------------------------

init.causalForest <- function(x, y, w, num.trees) {
  
  num.obs <- nrow(x)                                                            # Number of samples
  
  trees <- vector("list", num.trees)                                            # Initialize a vector
  inbag <- matrix(0, num.obs, num.trees)                                        # Initialize a matrix of zeros
  
  causalForest <- list(trees = trees,                                           # Define causalForest as a list
                       x = x, y = y, w = w, 
                       ntree = num.trees, 
                       inbag = inbag)
  
  class(causalForest) <- "causalForest"                                         # Define a class
  causalForest
  
}