# fusing output csvs of SDMs
library(tidyverse)

ver <- "VM3.1" # designate version

# list model output files
files <- list.files(paste0("/gpfs/gibbs/pi/jetz/from_loomis/data/results_30x30/birds/",
                           ver), full.names=T,  recursive=F)
files <- files[grep(".csv", files)]
# loop through each species file and combine
for (i in 1:length(files)){
  file <- read_csv(files[i])
  if(ncol(file)==64){file = cbind(file, "NA")} 
  colnames(file)[65] = "ghow"
  if(i>1){names(table) = names(file)}
  if(i==1){table <- file}else{table <- rbind(table, file)}
}
# write table summarizing all models
write_csv(table, paste0("/gpfs/gibbs/pi/jetz/from_loomis/data/results_30x30/birds/sdm_birds_",ver,".csv"))
