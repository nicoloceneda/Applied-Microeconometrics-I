# ------------------------------------------------------------------------------
# Propensity forest (procedure 2)
# ------------------------------------------------------------------------------

propensityForest <- function(X, Y, W, num.trees,                                # {X}, {Y}, {W}, {B}
                             sample.size = floor(length(Y) / 10),               # {s}
                             # mtry = ceiling(ncol(X)/3),
                             nodesize = 1) {                                    # {k}
  
  # ----------------------------------------------------------------------------
  # Define the variables
  # ----------------------------------------------------------------------------
  
  # Check for missing values
  if (any(is.na(X)) || any(is.na(Y)) || any(is.na(W))) {                        
    stop("There are missing values in the input.")
  }
  
  # Additional input variables
  num.obs <-nrow(X)                                                             # Number of samples {n}
  sample.size <- min(sample.size, floor(num.obs / 2))                           # Sample size {s}
  
  # ----------------------------------------------------------------------------
  # Build the forest
  # ----------------------------------------------------------------------------
  
  # Train a causal tree on (X,W) pairs with i in I, where outcome is treatment assignment
  # - Pick splits by maximizing classification accuracy for treatment assignments
  # - Goal is to have observations within a leaf that have same treatment propensity
  
  causalForest.honest <- init.causalForest(X, Y, W, num.trees)                  # Initialize the class [F2]
  
  print("Building trees ...")
  
  for (tree.index in 1:num.trees) {
    
    print(paste("Tree", as.character(tree.index)))
    
    # Define the data
    full.idx <- sample.int(num.obs, sample.size, replace = FALSE)               # Sample {I} of size {s} without replacement
    tree.DF = data.frame(X = X, W = W)                                          # Training data {(X,W)}
    
    #---------------------------------------------------------------------------
    # Propensity tree
    #---------------------------------------------------------------------------
    
    # Recursive partitioning of propensity tree using W with i in I - Breiman et. al (1984) 
    tree.propensity <- rpart(W ~ .,                                             
                             data = tree.DF[full.idx,], 
                             method = "class", 
                             minbucket = nodesize)
    
    class(tree.propensity) <- "causalTree"                                      # Define a class causalTree
    
    attr(tree.propensity, "ylevels") = NULL                                     # Hack to make predict.causalTree do the right thing
    
    #---------------------------------------------------------------------------
    # Honest tree
    #---------------------------------------------------------------------------
    
    # Refit causal tree to get honest tree [F3]
    tree.honest <- refit.causalTree(tree.propensity, 
                                    newx = tree.DF[full.idx,], 
                                    newy = Y[full.idx], 
                                    treatment = W[full.idx])
    
    causalForest.honest$trees[[tree.index]] <- tree.honest
    causalForest.honest$inbag[full.idx, tree.index] <- 1
  }
  
  return(causalForest.honest)
}