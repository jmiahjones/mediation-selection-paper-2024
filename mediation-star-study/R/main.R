# main.R
# Description: The main data analysis pipeline.

print("Beginning pipeline for outcome 1...")

source("./R/data-preproc.R", echo=TRUE)
source("./R/data-step2.R", echo=TRUE)

sl_m_file = "./cache/sl_m_sel_alldata.RData"
if(file.exists(sl_m_file)){
  load(sl_m_file)
} else {
  source("./R/SL-estimate-M.R", echo=TRUE)
}

outcome_name <- "g8mathco"
Y <- Y1
sel_df <- data.frame(full_dat_noY, y=Y)
sl_y_file = "./cache/sl_y1_sel_alldata.RData"
if(file.exists(sl_y_file)){
  load(sl_y_file)
} else {
  source("./R/SL-estimate-Y.R", echo=TRUE)
}

# source("./R/sl-funs.R")
sel_result_file = "./cache/sel_result_y1_alldata.RData"
source("./R/selection-results.R", echo=TRUE)

print("Pipeline for outcome 1 complete! Check for warnings...")

rm(list=objects())

print("Beginning pipeline for outcome 2...")

source("./R/data-preproc.R", echo=TRUE)
source("./R/data-step2.R", echo=TRUE)

sl_m_file = "./cache/sl_m_sel_alldata.RData"
load(sl_m_file)

outcome_name <- "g8math_a"
Y <- Y2
sel_df <- data.frame(full_dat_noY, y=Y)
sl_y_file = "./cache/sl_y2_sel_alldata.RData"
if(file.exists(sl_y_file)){
  load(sl_y_file)
} else {
  source("./R/SL-estimate-Y.R", echo=TRUE)
}

sel_result_file = "./cache/sel_result_y2_alldata.RData"
source("./R/selection-results.R", echo=TRUE)

print("Pipeline for outcome 2 complete! Check for warnings...")



# 
# rm(list=objects())
# 
# print("Beginning pipeline for outcome 3...")
# 
# source("./R/data-preproc.R", echo=TRUE)
# source("./R/data-step2.R", echo=TRUE)
# 
# sl_m_file = "./cache/sl_m_sel_alldata.RData"
# load(sl_m_file)
# 
# outcome_name <- "g4mathconcapplss"
# Y <- Y3.select
# sel_df <- data.frame(full_dat_noY, y=Y)
# sl_y_file = "./cache/sl_y3_sel_alldata.RData"
# if(file.exists(sl_y_file)){
#   load(sl_y_file)
# } else {
#   source("./R/SL-estimate-Y.R", echo=TRUE)
# }
# 
# sel_result_file = "./cache/sel_result_y3_alldata.RData"
# source("./R/selection-results.R", echo=TRUE)
# 
# print("Pipeline for outcome 2 complete! Check for warnings...")

