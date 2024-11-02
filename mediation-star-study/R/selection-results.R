## Inputs: outcome variable as Y
##         path to save file as sel_result_file

# label output by outcome
opt <- list(outcome=outcome_name)

# load functions for doing selection + estimation + inference
source("./R/mediation-funs-postsel.R", echo=FALSE)

# after running data through step 2 + SL-estimate-M + "Y choice" + SL-estimate-Y
mu.hats=c(
  mu_d_list,
  sl_m,
  sl_y
)

num_bootstraps = 2000L
lambda.len <- 401L
lambdas = n.dat^(-0.75)*(2^seq(from=-2, to=14, length.out=lambda.len))

weight.kappas <- c(0.5,1,2,3)

d <- pull(sel_df, treat.col)
m <- as.matrix(dplyr::select(sel_df, one_of(candidate.mediators)))
x <- as.matrix(dplyr::select(sel_df, one_of(x.cols)))
y <- pull(sel_df, outcome.col)

dc <- d - mu.hats$mu.dx
m_0 <- m - mu.hats$mu.mxis
y_0 <- y - mu.hats$mu.yx

print("Beginning prd")
prd_results <- cv_sl_estimates_w_sel(
  y_0, m_0, dc, weight.version="product", weight.gam=weight.kappas,
  lambdas=lambdas, folds=folds, opt=opt, 
  num_bootstraps=num_bootstraps, do.boot=TRUE, boot.seed=2021
)

print("Beginning adp")
adp_results <- cv_sl_estimates_w_sel(
  y_0, m_0, dc, weight.version="adaptive", weight.gam=weight.kappas,
  lambdas=lambdas, folds=folds, opt=opt, 
  num_bootstraps=num_bootstraps, do.boot=TRUE, boot.seed=2021
)

print("Beginning full")
full_results <- cv_sl_estimates_no_sel(
  y_0, m_0, dc, model_name="full", opt=opt, 
  num_bootstraps=num_bootstraps, do.boot=TRUE, boot.seed=2021
)


y_0_lm <- lm(y ~ x)$residuals
m_0_lm <- lm(m ~ x)$residuals
dc_lm <- lm(d ~ x)$residuals
lm_results <- cv_sl_estimates_no_sel(
  y_0_lm, m_0_lm, dc_lm,
  model_name="parametric",
  opt=opt,
  num_bootstraps=num_bootstraps,
  do.boot=TRUE,
  boot.seed=2021
)


save(prd_results, adp_results, full_results, lm_results,
     file=sel_result_file)