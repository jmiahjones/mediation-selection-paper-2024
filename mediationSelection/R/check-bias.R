#
# Author: Jeremiah Jones
# rm(list=objects())

check_bias <- function(
  n,
  num_simulations,
  abbv_scn,
  coef_setting,
  use_sl,
  suffix_arg,
  cores
) {
  
  assertthat::assert_that(is.logical(use_sl))
  
  ######################## Load Procedures ########################
  # need this in the function environment
  source("./R/simulation-helpers.R", local=T)
  source("./R/mediation-funs-postsel.R", local=T)
  
  num_noise = 7
  
  num_bootstraps=1000L
  small.setting <- tolower(coef_setting)=="small"
  small.alpha.setting <- tolower(coef_setting)=="smallalpha"
  many.p.setting <- tolower(coef_setting)=="manyp"
  is.randomized <- substr(abbv_scn, 1,1)=="r"
  
  if(small.setting){
    # small
    # chose alphas*betas=16/sqrt(n) for the first 3 mediators
    alphas = 4*c(n^(-1/4), 1, 1, rep(0, num_noise)) # D-M coef
    betas = c(n^(-1/4), n^(-1/2), n^(-1/2), rep(0,num_noise)) # M-Y coef
  } else if(small.alpha.setting) {
    alphas = n^(-1/2)*c(1, 1, n^(1/4), rep(0, num_noise)) # D-M coef
    betas = 16*c(1, 1, n^(-1/4), rep(0,num_noise)) # M-Y coef
  } else if(many.p.setting) {
    # large fixed coefficients
    num_noise <- 57 # p=60
    alphas <- c(1, 2, 2, rep(0, num_noise))
    betas <- c(.8, .4, .4, rep(0, num_noise))
  } else {
    # large fixed coefficients
    alphas <- c(1, 2, 2, rep(0, num_noise))
    betas <- c(.8, .4, .4, rep(0, num_noise))
  }
  
  
  
  noise = c(0,0,0,rep(0,num_noise))
  variances <- c(1,1,1, rep(1, num_noise-1), 1)
  NDE = 2
  NIE = sum(alphas * betas)
  rho = 0.0 # covariance between errors
  rho2 = 0.0 # noise covariance
  corr_meds = c(1,9)
  p = length(alphas)
  covariates_size = 3
  V = 10 # number of folds for crossvalidation
  
  
  
  ################# Initialization ########################
  
  start = Sys.time()
  
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(SuperLearner)
  library(glmnet)
  # library(gam)
  library(mgcv)
  library(earth)
  library(randomForest)
  # library(e1071)
  # library(xgboost)
  library(parallel)
  library(doParallel)
  library(foreach)
  # print(sessionInfo())
  
  
  candidate.mediators = paste0("m.", 1:p)
  true_M <- which(alphas*betas != 0)
  oracle.mediators <- candidate.mediators[true_M]
  m.cols="m"
  # x.cols = paste0("x.", 1:covariates_size)
  x.cols="x"
  expit <- plogis
  
  # Setup a lambda grid
  lambda.len = 301
  lambdas = n^(-0.75)*(2^seq(from=-2, to=10, length.out=lambda.len))
  
  set.seed(841664)
  folds <- caret::createFolds(y=1:n, k=V)
  
  ################# Files ########################
  
  savefile_suffix <- paste(n, num_simulations, abbv_scn, coef_setting,
                           use_sl, suffix_arg, sep="-")
  
  logfile <- paste0("./logs/cluster-", savefile_suffix, ".out")
  
  if(!dir.exists("./logs")){
    dir.create("./logs")
  }
  if(file.exists(logfile)){
    file.remove(logfile)
  }
  
  
  
  #################### Confounding Mechanisms ########################
  
  
  if(substr(abbv_scn, 1,1)=="n"){
    # nonlinear function
    propensity <- function(x) {
      x1 <- x[1]
      x2 <- x[2]
      x3 <- x[3]
      x4 <- x[4]
      x5 <- x[5]
      # expit((0*x1 + 1*x2 + 0*x3 + 0*x4 + 0*x5+0*x1^2*x2+ 1*x1*x3 - 0)/1)
      expit((x2 + x1*x2)/1.25)
    }
  } else if(substr(abbv_scn, 1,1)=="l"){
    # linear function
    propensity <- function(x) {
      x1 <- x[1]
      x2 <- x[2]
      x3 <- x[3]
      x4 <- x[4]
      x5 <- x[5]
      expit((1*x1 + 1*x2 )/1.25)
    }
  } else if(substr(abbv_scn, 1,1)=="r"){
    # linear function
    propensity <- function(x) {
      0.5
    }
  }
  
  
  if(substr(abbv_scn, 2,2)=="n"){
    # nonlinear function
    psi_m <- function(x) { 
      x1 <- x[1]
      x2 <- x[2]
      x3 <- x[3]
      x4 <- x[4]
      x5 <- x[5]
      1*x1^2 + x2 - x3 # 0*x3 - 0*x4 - 0*x5 + 0*x2*(x1-0.5)^2
    }
  } else if(substr(abbv_scn, 2,2)=="l"){
    # linear function
    psi_m <- function(x) { 
      x1 <- x[1]
      x2 <- x[2]
      x3 <- x[3]
      x4 <- x[4]
      x5 <- x[5]
      1*x1 + x2 - x3
    }
  }
  
  
  
  if(substr(abbv_scn, 3,3)=="n"){
    # nonlinear function
    psi_y <- function(x) {
      x1 <- x[1]
      x2 <- x[2]
      x3 <- x[3]
      x4 <- x[4]
      x5 <- x[5]
      1*(2*(x1-0.5)^2 + x2 + 2*x3)
    }
  } else if(substr(abbv_scn, 3,3)=="l"){
    # linear function
    psi_y <- function(x) {
      x1 <- x[1]
      x2 <- x[2]
      x3 <- x[3]
      x4 <- x[4]
      x5 <- x[5]
      (2*(x1-0.5) + x2 + 2*x3)
    }
  }
  
  
  ######################## Simulations ########################
  
  if(cores > 1){
    cl <- parallel::makeCluster(cores, "FORK", outfile=logfile)
    doParallel::registerDoParallel(cl)
    parallel::clusterSetRNGStream(cl, 2018)
  } else {
    foreach::registerDoSEQ()
  }
  
  
  
  # tryCatch({
  simulations <- foreach(sim.idx=1:num_simulations) %dopar% {
    
    #### Simulation mechanism ####
    
    # set seed to ensure similar xdm each time -- allows for saving sl
    set.seed(2021 + sim.idx)
    # confounders
    # x <- replicate(covariates_size, runif(n, min=0, max=1))
    x <- replicate(covariates_size, rnorm(n, sd=0.5))
    
    # treatment confounding
    d.prob <- apply(x, 1, function(row){
      propensity <- propensity(row)
    })
    d <- rbinom(n=n, size=1, prob=d.prob)
    
    Sigma <- 
      matrix(c(rep(rho, p**2)), ncol=p)
    
    Sigma[corr_meds, p] <- rho2
    Sigma[p, corr_meds] <- rho2
    Sigma[corr_meds, corr_meds] <- rho2
    diag(Sigma) <- rep(1, p)
    
    Sigma <- diag(sqrt(variances)) %*% Sigma %*%
      diag(sqrt(variances))
    # if(p>5){
    #   Sigma[,as.logical(noise)] <- rho2
    #   Sigma[as.logical(noise),] <- rho2
    # }
    
    epsilon_m <- MASS::mvrnorm(n, mu=rep(0, p), Sigma=Sigma)
    
    # mediator confounding
    confounding_m <- apply(x, 1, function(row){
      psi_m(row)
    })
    # browser()
    # mediators
    m <- foreach(i=1:n, .combine=rbind) %do%{
      d[i] * alphas + confounding_m[i] + epsilon_m[i,]
    }
    # m[,p] <- m[,-p] %>% rowMeans
    # m[,p] <- m[,1]^2 + m[,2]^2 + rowMeans(m[,4:(p-1)])
    # m <- m + epsilon_m
    
    # outcome confounding
    confounding_y <- apply(x, 1, function(row){
      psi_y(row)
    })
    
    # set.seed(7243957 + sim.idx)
    # outcome equation
    epsilon_y <- rnorm(n, mean=0, sd=1)
    y <- (NDE * d) + (m %*% betas) + (confounding_y) + epsilon_y
    
    # for troubleshooting only
    true_em <- foreach(i=1:n, .combine=rbind) %do%{
      (d.prob[i] * alphas) + confounding_m[i]
    }
    true_ey <- as.numeric(NDE*d.prob + true_em %*% as.matrix(betas) + confounding_y)
    
    stopifnot(max(abs(
      true_ey - confounding_y - NDE*d.prob - true_em %*% betas
    )) < 1e-10)
    
    # if(sim.idx %% 10 == 0)
    #   print(paste("Completed simulation", sim.idx))
    
    # save the true centering variables for future use
    # em = replicate(p,confounding_m) + replicate(p,propensity)*alphas
    # ey = NDE*propensity + (em %*% betas) + confounding_y
    
    # return(
    tibble(
      x=x, d=d, m=m, y=y,
      propensity=d.prob,
      true_em = true_em,
      true_ey = true_ey,
      psi_mx = confounding_m,
      psi_yx = confounding_y,
      epsilon = epsilon_y,
      epsilon_m = epsilon_m
    )
    # )
  }
  # return(simulations)
  print("Created data.")
  # parallel::clusterExport(cl, "simulations")
  
  
  opt <- tibble(
    # n = n, # already calculated
    num_simulations = num_simulations,
    scenario = abbv_scn,
    coef_setting = coef_setting,
    use_sl = use_sl,
    suffix = suffix_arg,
    cores = cores
  )
  
  #### Use the methods ####
  results <- foreach(sim.idx=1:num_simulations, sim=simulations,
                     .combine=rbind, .errorhandling="remove"
  ) %dopar% {
    
    fwl_d <- lm(sim$d ~ sim$x)
    fwl_m <- lm(sim$m[,true_M] ~ sim$x)
    fwl_y <- lm(sim$y ~ sim$x)
    
    dc <- fwl_d$residuals
    m_0 <- fwl_m$residuals
    y_0 <- fwl_y$residuals
    
    mu_err <- c(
      muderr=0,
      mumerr=rep(0, ncol(sim$m)),
      muyerr=0
    )
    
    names(mu_err)[-c(1, length(mu_err))] <- paste0("mumerr", 1:p)
    this_opt <- tibble(opt, mu_err %>% as.list %>% as_tibble)
    
    prd_results <- mix_results <- adp_results <- NULL
    
    print("Beginning parametric")
    oracle_results <- cv_sl_estimates_no_sel(
      y_0, m_0, dc, model_name="oracle-lm", opt=this_opt, 
      num_bootstraps=num_bootstraps, do.boot=TRUE, boot.seed=2349871
    ) %>% mutate(sim=sim.idx)
    
    return(
      oracle_results %>% mutate(
        coverage_NDE = (lower_NDE <= NDE) & (NDE <= upper_NDE),
        coverage_NIE = (lower_NIE <= NIE) & (NIE <= upper_NIE),
        err_NDE = NDE_hat - NDE,
        err_NIE = NIE_hat - NIE
      )
    )
  }
  
  
  print("Scenario Complete!")
  stop = Sys.time()
  print(stop - start)
  
  ret <- results
  return(ret)
  
}


library(qs)
library(foreach)

################# Call Main ########################
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  warning("No commandline arguments found. Debugging...")
  args <- c(
    100, # n
    20, # num_simulations
    "lll", # abbv_scn
    "large", # coef_setting
    "F", # use_sl
    "debug", # suffix_arg
    2 # cores
  )
}

DEFAULT_CORES <- 20

stopifnot(length(args) %in% (6:7))
n <- as.numeric(args[1])
num_simulations <- as.numeric(args[2])
abbv_scn <- args[3]
coef_setting <- args[4]
use_sl <- args[5]
suffix_arg <- args[6]
cores <- as.numeric(args[7])

assertthat::assert_that(use_sl %in% c("F","T"))
use_sl_str <- use_sl
use_sl <- use_sl == "T"

if(is.na(cores)){
  cores <- DEFAULT_CORES
}

savefile_suffix <- paste(n, num_simulations, abbv_scn, coef_setting,
                         use_sl_str, suffix_arg, sep="-")

save_file <- paste0("./results/result-",
                    savefile_suffix,
                    ".qs"
)

print(paste0("Cores: ", cores))

sl_loop <- F
main_result <- foreach(slB=sl_loop, .combine=rbind) %do%{
  check_bias(
    n,
    num_simulations,
    abbv_scn,
    coef_setting,
    slB,
    suffix_arg,
    cores
  )
}

# save the result in RDS
qs::qsave(main_result, save_file)


