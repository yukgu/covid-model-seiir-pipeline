###############################################################
# Make scatters
#              To run: use interactively and manually change as needed
#  1 - scatter a single run given two dates
#  2 - scatter two runs given a single date
###############################################################
rm(list = ls())
library(data.table)
library(ggplot2)
library(ggrepel)
library(viridis)
library(gridExtra)
source("/ihme/cc_resources/libraries/current/r/get_location_metadata.R")
source(paste0("/ihme/code/covid-19/user/",Sys.info()['user'],"/hospital-r/code/collapse_combine_functions.R"))
seir_root <- "/ihme/covid-19/hospitalizations/inputs/seir/"

# Set hierarchy and agg locations
hierarchy <- get_location_metadata(111, 771, release_id = 9)
agg_locations <- hierarchy[most_detailed==0 & level >= 2]
agg_locations <- agg_locations[order(level)]
#removing E Asia which has no data and China
agg_locations<-agg_locations[!location_id %in% c(5,6)]
compare_versions_plot <- function(version_1, version_2, agg_locations, hierarchy, plot_date){
  
  data_1 <-  fread(paste0("/ihme/covid-19-2/seir-outputs/",version_1,"/reference/output_summaries/cumulative_deaths.csv"))
  data_2 <-  fread(paste0("/ihme/covid-19-2/seir-outputs/",version_2,"/reference/output_summaries/cumulative_deaths.csv"))
  
  
  plotting_data <- merge(data_1[date == as.Date(plot_date)], data_2[date == as.Date(plot_date)], by = c("location_id", "observed", "date"))
  plotting_data <- merge(plotting_data, hierarchy[,c("location_id","super_region_name","region_name","ihme_loc_id","map_id", "level","parent_id","local_id","location_name")], by="location_id", all.x=T)
  
  #plotting_data <- plotting_data[space == "Cumulative deaths"]
  
  plotting_data[, diff := (mean.y - mean.x)]
  plotting_data[, pc:= (mean.y - mean.x)/mean.x]
  plotting_data$above_25<-0
  plotting_data[abs(pc) > 0.25, above_25 := 1]
  plotting_data$above_45<-0
  plotting_data[abs(pc) > 0.45, above_45 := 1]

  
  min <- min(plotting_data[location_id==1]$mean.x, plotting_data[location_id==1]$mean.y)
  max <- max(plotting_data[location_id==1]$mean.x, plotting_data[location_id==1]$mean.y)
  globe_plot <- ggplot(plotting_data[location_id==1], aes(x = mean.x , y = mean.y))+
    geom_point()+geom_abline(aes(intercept =0, slope = 1))+theme_bw()+ylim(ymin= min, ymax = max)+xlim(xmin= min, xmax = max)+
    labs(x = "Current public reference", y = "New reference", title = paste0("Cumulative deaths compare, Global (", plot_date, ")"))
  print(globe_plot)
  
  national_plot <- ggplot(plotting_data[level ==3], aes(x = mean.x , y = mean.y))+
    geom_point(aes(col=factor(above_45)))+geom_abline(aes(intercept =0, slope = 1))+theme_bw()+geom_text_repel(data = plotting_data[level ==3 & above_45==1], aes(label = location_name), col="darkorchid2")+
    scale_color_manual("Change from public reference", values=c("black", "darkorchid2"), labels=c("<45% difference", ">45% difference"))+
    labs(x = "Current public reference", y = "New reference", title = paste0("Cumulative deaths compare, all national locations (", plot_date, ")"))
  #print(national_plot)
  
  national_plotl <- ggplot(plotting_data[level ==3], aes(x = log10(mean.x) , y = log10(mean.y)))+
    geom_point(aes(col=factor(above_45)))+geom_abline(aes(intercept =0, slope = 1))+theme_bw()+geom_text_repel(data = plotting_data[level ==3 & above_45==1], aes(label = location_name), col="darkorchid2")+
    scale_color_manual("Change from public reference", values=c("black", "darkorchid2"), labels=c("<45% difference", ">45% difference"))+
    labs(x = paste0("Current public reference", " (log)"), y = paste0("New reference", " (log)"), title = paste0("Cumulative deaths compare, all national locations (", plot_date, ")"))
  #print(national_plot)
  
  grid.arrange(national_plot,national_plotl, ncol=1)

  for(agg_loc in agg_locations$location_id){
    
    agg_name <- agg_locations[location_id == agg_loc, location_name]
    
    national_plot <- ggplot(plotting_data[parent_id == agg_loc], aes(x = mean.x , y = mean.y))+
      geom_point(aes(col=factor(above_25)))+geom_abline(aes(intercept =0, slope = 1))+theme_bw()+geom_text_repel(data=plotting_data[above_25==1 & parent_id==agg_loc],aes(label = location_name), col="cyan3")+
      scale_color_manual("Change from public reference", values=c("black", "cyan3"), labels=c("<25% difference", ">25% difference"))+
      labs(x = "Current public reference", y = "New reference", title = paste0("Cumulative deaths compare ",agg_name, " (", plot_date, ")"))
    #print(national_plot)
    
    national_plotl <- ggplot(plotting_data[parent_id == agg_loc], aes(x = log10(mean.x) , y = log10(mean.y)))+
      geom_point(aes(col=factor(above_25)))+geom_abline(aes(intercept =0, slope = 1))+theme_bw()+geom_text_repel(data=plotting_data[above_25==1 & parent_id==agg_loc],aes(label = location_name), col="cyan3")+
      scale_color_manual("Change from public reference", values=c("black", "cyan3"), labels=c("<25% difference", ">25% difference"))+
      labs(x = paste0("Current public reference", " (log)"), y = paste0("New reference", " (log)"), title = paste0("Cumulative deaths compare ",agg_name," (", plot_date, ")"))
    #print(national_plot)
    grid.arrange(national_plot,national_plotl, ncol=1)
  }
  
}
# Plot 2 compare two versions for a given date
version_1<-"2021_06_09.01" #Always the current public ref
version_2<-"2021_06_15.01" #new version

# Set plot_out
# plot_out <- paste0("/snfs1/Project/covid/results/diagnostics/seir/", version_2) #update to your desired location
# dir.create(plot_out)
# Sys.sleep(18000)
#dir.create(paste0(plot_out,version_2))
#Per Emm - only want last day of data for public viz - right now 2021-02-01
plot_out <- "/homes/aucarter/Downloads"
pdf(paste0(plot_out,"/seir_scatter_scenarios.pdf"), height=8, width=12)

# plot_date <- "2021-01-01"
# compare_versions_plot(version_1, version_2, agg_locations = agg_locations, hierarchy = hierarchy, plot_date = plot_date)

plot_date <- "2021-09-01"
compare_versions_plot(version_1, version_2, agg_locations = agg_locations, hierarchy = hierarchy, plot_date = plot_date)

# plot_date <- "2021-12-31"

dev.off()