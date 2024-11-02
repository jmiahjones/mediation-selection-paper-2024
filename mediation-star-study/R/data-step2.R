# data-step2.R
# Description: Create the selection set and inference set 
#              for subsequent analysis.

# Helpful to scale down the training set for testing.
selection_set_size <- 1

full_dat_noY <- data.frame(Xmat,
                           d=d_bin,
                           M_df)
x.cols <- colnames(Xmat)
treat.col <- "d"
candidate.mediators <- colnames(M_df)
outcome.col <- "y"

# partition selection set and inference set
n.dat <- nrow(full_dat_noY)


#### Estimate DX ####
d_sel <- dplyr::pull(full_dat_noY, treat.col)
mu_d <- mean(d_bin)
mu_d_list <- list(mu.dx=mu_d)

set.seed(2021)
folds <- caret::createFolds(d_sel, k=10L)
