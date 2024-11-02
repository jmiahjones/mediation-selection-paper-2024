library(glmnet)
library(tidyverse)
library(future.apply)
plan("multicore", workers=12L)
# plan(sequential)

start <- Sys.time()
message(paste0("Start time: ", start))

set.seed(2022, kind="L'Ecuyer-CMRG")

library(qs)
library(RcppArmadillo)
Rcpp::sourceCpp("coord_desc_arma.cpp")

.do_product_asym_cd <- function(d, M, Y, alpha_til, beta_til, lam_type, kappa){

  n <- nrow(M)

  # fixed-asymptotics rate
  if(lam_type == "fixed"){
    lam_star <- n**(1/6) # works for kappa >= 1
  }
  # local-asymptotics rate
  if(lam_type == "local"){
    lam_star <- n**(-(kappa-1)/2 - kappa/4)
  }

  pen_fac = abs(alpha_til * beta_til)**(-kappa)

  cd_coefs <- coord_desc(cbind(d,M), Y, w=c(0,pen_fac), lam=lam_star)

  list(best_coefs = cd_coefs |> as.numeric(),
       lambda = lam_star,
       kappa=kappa)
}

do_product_asym_cd <- function(d, M, Y, alpha_til, beta_til, lam_type, kappa, cv_type){
  res <- .do_product_asym_cd(d, M, Y, alpha_til, beta_til, lam_type, kappa)
  return(res)
}

do_sim_paths <- function(n, p, rho, use_lm, alpha_rate=0, beta_rate=0,
                         kappa=c(1,3), lam_type=c("local"),
                         cv_type=c(NA)) {

  stopifnot(all(kappa >= 1 | is.na(kappa))) # required based on asymp code
  stopifnot(all(lam_type %in% c("fixed", "local")))
  stopifnot(all(cv_type %in% c(NA,"mse", "bic")))
  stopifnot(alpha_rate <= 0)
  stopifnot(beta_rate <= 0)

  d <- rbinom(n, 1, .5)

  sigma11 <- matrix(rho, nrow=2, ncol=2)
  diag(sigma11) <- 1

  sigma22 <- matrix(rho, nrow=p-2, ncol=p-2)
  diag(sigma22) <- 1

  sigma <- matrix(0, nrow=p, ncol=p)
  sigma[1:2,1:2] <- sigma11
  sigma[3:p,3:p] <- sigma22
  sigma[3:4,1:2] <- sigma[1:2,3:4] <- 0.0

  eta <- MASS::mvrnorm(n=n, mu=rep(0,p), Sigma = sigma)

  M <- matrix(nrow=n, ncol=p)

  M[,1] <- n**(alpha_rate)*d + eta[,1]
  M[,2] <- d + eta[,2]
  M[,3] <- eta[,3]      # outcome predictor
  M[,4] <- d + eta[,4]  # effected by treatment
  M[,5:p] <- eta[,5:p]

  muy <- 2*d + n**(beta_rate) * M[,1] + rowSums(M[,2:3])

  Y <- muy + rnorm(n)

  alpha_til <- coef(lm(M ~ d))[2,]
  beta_til <- if(use_lm){
    coef(lm(Y ~ d + M))[-c(1:2)]
  } else {
    fit <- MASS::lm.ridge(Y~ d + M, lambda=2**seq(-10,4,length.out=501))
    lam <- which.min(fit$GCV)
    coef(fit)[lam,-c(1:2)]
  }

  configs <- expand_grid(kappa=kappa, lam_type=lam_type,
                         cv_type=cv_type) %>%
    filter((!is.na(kappa) & is.na(cv_type)) |
             (is.na(kappa) & !is.na(cv_type)))
  results <- configs %>% mutate(
    data = pmap(., \(...) do_product_asym_cd(..., d=d, M=M, Y=Y, alpha_til=alpha_til, beta_til=beta_til))
  )

  results$out <- lapply(seq.int(nrow(results)), \(res_idx) {
    res <- results$data[[res_idx]]
    best_coefs <- res$best_coefs
    # best_coefs <- best_coefs[-1]
    beta_hats <- best_coefs[-1]
    NIE_hat <- sum(alpha_til * beta_hats)
    NDE_hat <- best_coefs[1]

    Shat_B <- abs(beta_hats)>0

    return(
      tibble(
        ests=list(tibble(NIE_hat = NIE_hat,
                         NDE_hat = NDE_hat)),
        sel_result = list(tibble(
          Shat_B = Shat_B,
          sel_idx = seq_along(Shat_B)
        )),
        actual_kappa = res$kappa
      )
    )
  })
  results %>% select(-data)
}

simulate <- function(num_simulations, n, p, rho, use_lm, alpha_rate, beta_rate) {
  # tmp <- future_replicate(num_simulations, do_sim_paths(n, p, rho, use_lm, lam_type, kappa), simplify=F)
  tmp <- future_replicate(num_simulations, do_sim_paths(n, p, rho, use_lm, alpha_rate, beta_rate), simplify=F)
  tmp2 <- do.call(rbind, tmp) %>% unnest(out)
  # browser()
  configs <- c("kappa", "lam_type", "cv_type")

  sel_result <- tmp2 %>% select(kappa, lam_type, cv_type, sel_result) %>%
    unnest(sel_result) %>%
    mutate(sel_idx=if_else(sel_idx < 5, sel_idx, 5L)) %>%
    group_by(kappa, lam_type, cv_type, sel_idx) %>% summarize(Shat = mean(Shat_B)) %>% ungroup %>%
    group_by(kappa, lam_type, cv_type) %>% nest(.key="sel_result")

  ests <- tmp2 %>% select(kappa, lam_type, cv_type, ests) %>% unnest(ests) %>%
    group_by(kappa, lam_type, cv_type) %>%
    summarize(across(c(NDE_hat, NIE_hat),
                     list(mean=~mean(.x), sd=~sd(.x)))) %>% ungroup %>%
    group_by(kappa, lam_type, cv_type) %>% nest(.key="ests")

  kappa_summary <- tmp2 %>% select(kappa, lam_type, cv_type, actual_kappa) %>%
    count(kappa, lam_type, cv_type, actual_kappa) %>%
    ungroup %>% mutate(prob=n/num_simulations) %>%
    select(kappa, lam_type, cv_type, actual_kappa, prob) %>%
    group_by(kappa, lam_type, cv_type) %>% nest(.key="kappa_summary")

  message(sprintf("Finished simulation: %s", paste(num_simulations, n, p, rho, use_lm)))

  return(
    sel_result %>% inner_join(ests, by=configs) %>%
      inner_join(kappa_summary, by=configs)
  )
}

sim_params <- tribble(
  ~n,   ~p, ~rho, ~use_lm, ~alpha_rate, ~beta_rate,
  200,  100, .5, F, -1/4, -3/2,
  400,  100, .5, F, -1/4, -3/2,
)

full_results <- sim_params %>%
  mutate(
    data = pmap(., simulate, num_simulations=1000L)
  )

qsave(full_results, file = "./results/too_small_results.qs")

stop <- Sys.time()
message(paste0("Stop time: ", stop))
message(sprintf("Elapsed: %.2f minutes.", difftime(stop,start,units="mins")))

