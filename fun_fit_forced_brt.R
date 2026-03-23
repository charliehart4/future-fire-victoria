# Function for fitting BRTs with FORCED parameters

# Function modified version of function by Charlie Hart for RAC
fit_forced_brt = function(data.in, initialLearningRate, seed, max_trees = 20000, treeComplexity = 2, bagFraction = 0.6){
  
  # Initialize values
  ntrees <- 0
  modelBRT <- NULL
  counter <- 0
  learningRate <- initialLearningRate
  
  # Fit model with set seed using forced parameters
  message(paste("Fitting BRT with Forced Parameters | LR:", learningRate, "| TC:", treeComplexity, "| BF:", bagFraction, "| Max Trees:", max_trees))
  
  # Fit model with set seed
  withr::with_seed(seed, {
    modelBRT <-
      dismo::gbm.step(data = data.in,
                      gbm.x = 4:ncol(data.in),
                      gbm.y = 3,
                      family = "bernoulli",
                      tree.complexity = treeComplexity,
                      #fold.vector = data.in$spatial_temporal_folds,
                      n.folds = 15,
                      #length(unique(data.in$spatial_temporal_folds)),
                      learning.rate = learningRate,
                      bag.fraction = bagFraction,
                      max.trees = max_trees,
                      plot.main = FALSE)})
  
  
  if (is.null(modelBRT)) {
    warning("Failed to fit a valid model after 1 attempt.")
  } else {
    message("Final model fitted with learning rate:", learningRate)
    return(modelBRT)
  }
}