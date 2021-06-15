# PURPOSE: create a csv file of cumulative deaths for each location, comparing
# test run (all scenarios) to current public reference
#
# INSTRUCTIONS: Specify SEIR version + scenario name to indicate which scenarios you want to include in the csv.
# You should also specify which dates to include. The dates should always be: (1) data date (2) end date for public
# forecasts (3) one month prior to the current end date (4) Dec 31, 2021


##### SET UP
rm(list = ls())
library(data.table)
library(yaml)
library(dplyr)
source("/ihme/cc_resources/libraries/current/r/get_location_metadata.R")

data_dir <- '/ihme/covid-19/seir-forecast/'
data_dir2 <- '/ihme/covid-19/seir-outputs/'

who <- 3; euro <- 118; wb <- 117; gbd <- 111

##### SET ARGS

# Which SEIR version(s) do you want to include in the csv?
seir_versions <- c('2021_06_09.01',
                   '2021_06_15.01'
)
# Which SEIR scenarios (from the above SEIR versions) do you want to include in the csv?
# We don't expect this to change
scenario_names <- c('reference',
                    'reference')

#How do you want to label the columns within the csv?
#We don't expect this to change
scenario_labels <-  c('public_ref',
                      'gbd_ref')

# Which dates do you want to include?
dates <- c('2021-06-01', # data date - update this to pull automatically from metadata.yaml file
           '2021-09-01', #tentative end date
           '2021-12-31'  # end date for internal reporting
)

# Where do you want to output the csv?
out_date<-as.Date(dates[1],format="%Y-%m-%d")
# save_dir <- paste0("/snfs1/Project/covid/results/diagnostics/seir/",seir_versions[5])
# dir.create(save_dir)
# system(paste0("chmod -R g+rwx ", save_dir))
save_dir <- "/homes/aucarter/Downloads/"
### Build the CSV
dates <- as.Date(dates)
date_labels <- format(dates, "%b%d_%Y")

# (1) SLT vetting csv ----------------------------------------------------------------------------------------------------------

all_data <- data.table()
i <- 1
for(v in seir_versions){
  scenario <- scenario_names[i]
  print(paste(v,scenario,"scenario"))
  label <- scenario_labels[i]
  print(label)
  data <- fread(paste0(data_dir2,v,'/',scenario,'/output_summaries/cumulative_deaths.csv'))
  data2<-fread(paste0(data_dir2,v,'/',scenario,'/output_summaries/cumulative_unscaled_deaths.csv'))
  setnames(data2, 'mean','unscaled')
  data<-setDT(merge(data, data2[,.(location_id, date, unscaled)], by=c('location_id','date')))
  for (d in 1:length(dates)) {
    tmp_deaths <- data[date == dates[d],][, paste0(date_labels[d],"_", label):= mean]
    tmp_deaths <- tmp_deaths[date == dates[d],][, paste0(date_labels[d],"_", label,'_unscaled'):=unscaled]
        if(d==1){
      #if(i!=1){tmp_data <- NULL}
      tmp_data <- tmp_deaths[,c(1,8,9)]
    }else{
      tmp_data <- merge(tmp_data, tmp_deaths[,c(1,8,9)], by = "location_id")
    }
  }

  if(i==1){
    all_data <- tmp_data
  }else{
    all_data <- merge(all_data, tmp_data, by = "location_id", all = T)
  }
  i <- i + 1
}

# Calculate diff and pct diff (current public ref vs ref)
all_data <- all_data[, diff_prev_current := get(paste0(date_labels[1],"_gbd_ref")) - get(paste0(date_labels[1],"_public_ref"))]
all_data <- all_data[, pct_diff_prev_current := (get(paste0(date_labels[1],"_gbd_ref")) - get(paste0(date_labels[1],"_public_ref")))/get(paste0(date_labels[1],"_public_ref"))]


# output csv
# write.csv(alldata, paste0(save_dir, "/cumulative_death_compare_seir_",seir_versions[5],".csv"), row.names = F)
write.csv(all_data, "/homes/aucarter/Downloads/cumulative_death_compare_seir.csv", row.names = F)
