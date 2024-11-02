# define-SL.R
library(SuperLearner)
library(nnet)
library(randomForest)
library(xgboost)
library(gam)
library(earth)
library(glmnet)

# nnet_lrn <- SuperLearner::create.Learner("SL.nnet",
#                                          tune=list(size=c(5,10,15,20))
# )
# nnet_names <- nnet_lrn$names

rf_lrn <- SuperLearner::create.Learner("SL.randomForest",
                                       params=list(ntree=1000,
                                                   mtry=3),
                                       tune=list(nodesize=c(5,9))
)
rf_names <- rf_lrn$names

SL.xgboost <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000, 
          max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(), 
          nthread = 1, verbose = 0, save_period = NULL, ...) 
{
  
  if (!is.matrix(X)) {
    X = model.matrix(~. - 1, X)
  }
  xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
  if (family$family == "gaussian") {
    model = xgboost::xgboost(data = xgmat, objective = "reg:squarederror", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period)
  }
  if (family$family == "binomial") {
    model = xgboost::xgboost(data = xgmat, objective = "binary:logistic", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period)
  }
  if (family$family == "multinomial") {
    model = xgboost::xgboost(data = xgmat, objective = "multi:softmax", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, num_class = length(unique(Y)), 
                             nthread = nthread, params = params, save_period = save_period)
  }
  if (!is.matrix(newX)) {
    newX = model.matrix(~. - 1, newX)
  }
  pred = predict(model, newdata = newX)
  fit = list(object = model)
  class(fit) = c("SL.xgboost")
  out = list(pred = pred, fit = fit)
  return(out)
}

xgb_lrn <- SuperLearner::create.Learner("SL.xgboost",
                                        params = list(max_depth=c(3)),
                                        tune=list(ntrees=c(100))
)
xgb_names <- xgb_lrn$names

cont_lib <- c(
  "SL.lm",
  "SL.gam",
  "SL.glmnet",
  "SL.earth",
  rf_names,
  xgb_names,
  # nnet_names,
  "SL.mean"
)
