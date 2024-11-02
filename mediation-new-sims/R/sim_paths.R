library(glmnet)
library(tidyverse)
library(future.apply)
plan("multicore", workers=4L)
# plan(sequential)

set.seed(2022, kind="L'Ecuyer-CMRG")

library(qs)

do_product <- function(d, M, Y, alpha_til, beta_til){
  S_star_B <- c(F,T,T,T, rep(F,ncol(M)-4))
  lambda_oracle <- NA
  oracle_coefs <- NULL
  best_fit <- NULL
  best_kap <- NULL
  min_risk <- Inf
  kappas <- c(.5,1,2,3,4,5)
  folds <- as.integer(cut(1:length(Y), breaks=10L))
  lambdas <- 2**seq(-10, 6, length.out=501)[-1]
  for(kappa in kappas) {
    pen_fac = abs(alpha_til * beta_til)**(-kappa)
    cvfit=cv.glmnet(cbind(d,M), Y,
                    penalty.factor = c(0, pen_fac), foldid = folds,
                    lambda=lambdas)
    this_risk <- min(cvfit$cvm)
    if(this_risk < min_risk){
      min_risk <- this_risk
      best_fit <- cvfit
      best_kap <- kappa
    }

  }

  list(
    best_coefs = coef(best_fit, s="lambda.min") |> as.numeric(),
    oracle_coefs = oracle_coefs
  )
}

do_sim_paths <- function(n, p, m1_size, use_lm) {

  d <- rbinom(n, 1, .5)

  eta <- replicate(p, rnorm(n))

  M <- matrix(nrow=n, ncol=p)

  M[,1] <- m1_size*d + eta[,1]
  M[,2] <- M[,1] + eta[,2]
  M[,3] <- m1_size*d + eta[,3]
  M[,4] <- M[,3] + eta[,4]
  M[,5:p] <- eta[,5:p]

  muy <- 2*d + .8*M[,2] + M[,3] + M[,4]

  Y <- muy + rnorm(n)

  S_star_B <- c(F,T,T,T, rep(F,p-4))

  alpha_til <- coef(lm(M ~ d))[2,]
  beta_til <- if(use_lm){
    coef(lm(Y ~ d + M))[-c(1:2)]
  } else {
    fit <- MASS::lm.ridge(Y~ d + M, lambda=2**seq(-10,4,length.out=501))
    lam <- which.min(fit$GCV)
    coef(fit)[lam,-c(1:2)]
  }

  conv_idx <- 0L
  not_converged <- T
  res <- do_product(d, M, Y, alpha_til, beta_til)
  while(not_converged && conv_idx < 5L){
    conv_idx <- conv_idx + 1L
    last_res <- res
    res <- do_product(d, M, Y, alpha_til, res$best_coefs[-(1:2)])
    not_converged <- mean(abs(res$best_coefs - last_res$best_coefs)>1e-8)
  }
  # if(not_converged) warning("Terminated without convergence.")
  # message("Terminated after ", conv_idx, " iterations.")


  # message(sprintf("Best kappa=%.1f, Best lambda=%.2f", best_kap,
  #         best_fit$lambda.min))

  best_coefs <- res$best_coefs
  best_coefs <- best_coefs[-1]
  beta_hats <- best_coefs[-1]
  NIE_hat <- sum(alpha_til * beta_hats)
  NDE_hat <- best_coefs[1]

  Shat_B <- abs(beta_hats)>0

  return(
    tibble(ests=list(tibble(NIE_hat = NIE_hat,
                NDE_hat = NDE_hat)),
           sel_result = list(tibble(
             Shat_B = Shat_B,
             sel_idx = seq_along(Shat_B)
           ))
    )
  )
}

simulate <- function(num_simulations, n, p, m1_size, use_lm) {
  tmp <- future_replicate(num_simulations, do_sim_paths(n, p, m1_size, use_lm), simplify=F)
  sel_result <- do.call(rbind, tmp) %>% select(-ests) %>% unnest(sel_result) %>%
    mutate(sel_idx=if_else(sel_idx < 5, sel_idx, 5L)) %>%
    group_by(sel_idx) %>% summarize(Shat = mean(Shat_B)) %>% ungroup

  ests <- do.call(rbind, tmp) %>% select(ests) %>% unnest(ests) %>%
    summarize(across(everything(),
                     list(mean=~mean(.x), sd=~sd(.x))))

  message(sprintf("Finished simulation: %s", paste(num_simulations, n, p, m1_size, use_lm)))

  return(tibble(ests = list(ests), sel_result=list(sel_result)))
}

sim_params <- tribble(
  ~n,   ~p, ~m1_size, ~use_lm,#  ~num_simulations,
  200, 100, .5, F,
  200, 100, 1, F,
  200, 200, .5, F,
  200, 200, 1, F,

  200, 100, 2, F,
  200, 100, 4, F,
  200, 200, 2, F,
  200, 200, 4, F,
)

full_results <- sim_params %>%
  mutate(
    data = pmap(., simulate, num_simulations=1000L)
  )


qsave(full_results, file = "./results/full_results.qs")
