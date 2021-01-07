# ------------------------------------------------------------------------------
# Initialize propensity forest
# ------------------------------------------------------------------------------

init.causalForest <- function(x, y, w, num.trees) {
  
  # Define the variables to include in the causalForest
  trees <- vector("list", num.trees)                                            # Initialize a vector of trees
  num.obs <- nrow(x)                                                            # Number of samples {n}
  inbag <- matrix(0, num.obs, num.trees)                                        # Initialize a matrix of zeros
  
  # Define causalForest as a list
  causalForest <- list(trees = trees,                                           
                       x = x, y = y, w = w, ntree = num.trees, 
                       inbag = inbag)
  # Create a class
  class(causalForest) <- "causalForest"                                         
  causalForest
}