## SL-estimate-Y.R
## Inputs: outcome variable as Y
##         path to save file as sl_y_file

source("./R/mediation-funs-postsel.R", echo=FALSE)
source("./R/sl-funs.R")
source("./R/define-SL.R")
# returns cont_lib

sl_y <- train_sl_mu(sel_df, m.cols=candidate.mediators, 
                    x.cols=x.cols, treat.col=treat.col, 
                    outcome.col=outcome.col,
                    bin_lib=NULL, 
                    cont_lib=cont_lib, 
                    folds=folds, 
                    cores=1L,
                    seed=2021L,
                    parallel_outfile="./logs/cluster-y.log",
                    save_pred_only=FALSE,
                    estD=F, 
                    estM=F,
                    estY=T,
                    verbose=T)

save(sl_y, file=sl_y_file)

