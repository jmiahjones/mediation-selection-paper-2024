## SL-estimate-M.R
# Inputs: save directory as sl_m_file

source("./R/mediation-funs-postsel.R", echo=FALSE)
source("./R/sl-funs.R")
source("./R/define-SL.R")
# returns cont_lib

sl_m <- train_sl_mu(full_dat_noY, m.cols=candidate.mediators, 
                    x.cols=x.cols, treat.col=treat.col, 
                    outcome.col=outcome.col,
                    bin_lib=NULL, 
                    cont_lib=cont_lib, 
                    folds=folds, 
                    cores=5L,
                    seed=2021L,
                    parallel_outfile="./logs/cluster-m.log",
                    save_pred_only=FALSE,
                    estD=F, 
                    estM=T,
                    estY=F,
                    verbose=T)

save(sl_m, file=sl_m_file)

