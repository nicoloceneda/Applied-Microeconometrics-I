# ------------------------------------------------------------------------------
# Propensity Forest
# ------------------------------------------------------------------------------

# Define propensity trees (procedure 2)
# - Use only the treatment indicator W to place splits
# - Use the responses Y for estimating tau

propensityForest <- function(X, Y, W, 
                             num.trees,                                         # Number of trees {B}
                             sample.size = floor(length(Y) / 10),               # Sample size {s}
                             mtry = ceiling(ncol(X)/3),
                             nodesize = 1) {                                    # Minimum leaf size {k}
  
  # Check for missing values
  if (any(is.na(X)) || any(is.na(Y)) || any(is.na(W))) {                        
    stop("There are missing values in the input.")
  }
  
  # Define the variables
  num.obs <-nrow(X)                                                             # Number of samples {n}
  causalForest.honest <- init.causalForest(X, Y, W, num.trees)                  # Initialize the class
  sample.size <- min(sample.size, floor(num.obs / 2))                           # Sample size {s}
  
  print("Building trees ...")
  
  # Train a causal tree on (X,W) pairs with i in I, where outcome is treatment assignment
  # - Pick splits by maximizing classification accuracy for treatment assignments
  # - Goal is to have observations within a leaf that have same treatment propensity
  
  for (tree.index in 1:num.trees) {
    
    print(paste("Tree", as.character(tree.index)))
    
    # Define the data
    full.idx <- sample.int(num.obs, sample.size, replace = FALSE)               # Sample {I} of size {s} without replacement
    tree.DF = data.frame(X = X, W = W)                                          # Training data {(X,W)}
    
    # Recursive partitioning of propensity tree using W with i in I
    tree.propensity <- rpart(W ~ .,                                             
                             data = tree.DF[full.idx,], 
                             method = "class", 
                             minbucket = nodesize)
    
    class(tree.propensity) <- "causalTree"                                      # Define a class
    
    attr(tree.propensity, "ylevels") = NULL                                     # Hack to make predict.causalTree do the right thing
    
    tree.honest <- refit.causalTree(tree.propensity, 
                                    newx = tree.DF[full.idx,], 
                                    newy = Y[full.idx], 
                                    treatment = W[full.idx])
    
    causalForest.honest$trees[[tree.index]] <- tree.honest
    causalForest.honest$inbag[full.idx, tree.index] <- 1
  }
  
  return(causalForest.honest)
}