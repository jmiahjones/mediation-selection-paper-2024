#!/bin/bash
# Name: main.sh
# Description: Starts the analysis pipeline in main.R by using a no-hangup
#              (nohup) call to Rscript, ensuring a clean session.

# uncomment to purge the cache
rm ./cache/*.RData

# uncomment to purge logs
rm ./logs/*.log

nohup Rscript ./R/main.R > ./logs/main.log
exit 0