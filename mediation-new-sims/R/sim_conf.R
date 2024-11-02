library(glmnet)
library(tidyverse)
library(future.apply)
plan("multicore", workers=10L)
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
  if(is.na(kappa)){
    kappa_grid <- c(1,3)
    nfolds <- 10L
    folds <- cut(1L:nrow(M), breaks=nfolds) |> factor() |> as.integer()

    tmp <- lapply(kappa_grid, \(kappa) {
        .do_product_asym_cd(
          d, M, Y,
          alpha_til, beta_til, lam_type, kappa)
    })

    if(cv_type=="mse"){
      emp_risks <- sapply(kappa_grid, \(kappa) {
        errs_sq_list <- lapply(seq.int(nfolds), \(foldid){
          foldidx <- folds == foldid
          cf <- .do_product_asym_cd(
            d[!foldidx], M[!foldidx,], Y[!foldidx],
            alpha_til, beta_til, lam_type, kappa)$best_coefs
          errs_sq <- (Y[foldidx] - cbind(d,M)[foldidx,] %*% cf)**2
          return(errs_sq)
        })
        mean(do.call(c, errs_sq_list))
      })
    } else {
      emp_risks <- sapply(kappa_grid, \(kappa) {
        loss_list <- lapply(seq.int(nfolds), \(foldid){
          foldidx <- folds == foldid
          cf <- .do_product_asym_cd(
            d[!foldidx], M[!foldidx,], Y[!foldidx],
            alpha_til, beta_til, lam_type, kappa)$best_coefs
          errs_sq <- (Y[foldidx] - cbind(d,M)[foldidx,] %*% cf)**2
          bic_corr <- log(nrow(M))*sum(abs(cf) > 1e-6)
          return(errs_sq + bic_corr)
        })
        mean(do.call(c, loss_list))
      })
    }
    min_idx <- which.min(emp_risks)
    res <- tmp[[min_idx]]
  } else {
    res <- .do_product_asym_cd(d, M, Y, alpha_til, beta_til, lam_type, kappa)
  }
  return(res)
}

do_sim_paths <- function(n, p, rho, use_lm, kappa=c(1,3,NA), lam_type=c("fixed","local"),
                         cv_type=c(NA,"mse","bic")) {

  stopifnot(all(kappa >= 1 | is.na(kappa))) # required based on asymp code
  stopifnot(all(lam_type %in% c("fixed", "local")))
  stopifnot(all(cv_type %in% c(NA,"mse", "bic")))

  d <- rbinom(n, 1, .5)

  # eta <- replicate(p, rnorm(n))
  sigma11 <- matrix(rho, nrow=2, ncol=2)
  diag(sigma11) <- 1

  sigma22 <- matrix(rho, nrow=p-2, ncol=p-2)
  diag(sigma22) <- 1

  sigma <- matrix(0, nrow=p, ncol=p)
  sigma[1:2,1:2] <- sigma11
  sigma[3:p,3:p] <- sigma22
  sigma[3:7,1:2] <- sigma[1:2,3:7] <- 0.25

  eta <- MASS::mvrnorm(n=n, mu=rep(0,p), Sigma = sigma)

  M <- matrix(nrow=n, ncol=p)

  M[,1] <- d + eta[,1]
  M[,2] <- d + eta[,2]
  M[,3:5] <- eta[,3:5]      # outcome predictor
  M[,6] <- d + eta[,6]  # effected by treatment
  M[,7:p] <- eta[,7:p]

  muy <- 2*d + rowSums(M[,1:5])

  Y <- muy + rnorm(n)

  alpha_til <- coef(lm(M ~ d))[2,]
  beta_til <- if(use_lm){
    coef(lm(Y ~ d + M))[-c(1:2)]
  } else {
    fit <- MASS::lm.ridge(Y~ d + M, lambda=2**seq(-10,4,length.out=501))
    lam <- which.min(fit$GCV)
    coef(fit)[lam,-c(1:2)]
  }

  # browser()
  # conv_idx <- 1L
  # conv_iter <- 1L
  # not_converged <- T

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

simulate <- function(num_simulations, n, p, rho, use_lm) {
  tmp <- future_replicate(num_simulations, do_sim_paths(n, p, rho, use_lm), simplify=F)
  tmp2 <- do.call(rbind, tmp) %>% unnest(out)
  configs <- c("kappa", "lam_type", "cv_type")

  sel_result <- tmp2 %>% select(kappa, lam_type, cv_type, sel_result) %>%
    unnest(sel_result) %>%
    mutate(sel_idx=if_else(sel_idx < 7, sel_idx, 7L)) %>%
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
  ~n,   ~p, ~rho, ~use_lm,
  200,  100, .5, F,
  400,  100, .5, F,
)

full_results <- sim_params %>%
  mutate(
    data = pmap(., simulate, num_simulations=1000L)
  )

qsave(full_results, file = "./results/conf_results.qs")

stop <- Sys.time()
message(paste0("Stop time: ", stop))
message(sprintf("Elapsed: %.2f minutes.", difftime(stop,start,units="mins")))

