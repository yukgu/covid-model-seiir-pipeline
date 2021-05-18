library(data.table); library(ggplot2); library(parallel)
root_dir <- "/ihme/covid-19-2/seir-outputs/"

# Parameters
output_version <- "2021_05_13.03"
quants <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
start_date <- as.Date("2021-05-02")
end_date <- as.Date("2022-03-01")

# Tables
source(file.path("/ihme/cc_resources/libraries/current/r/get_location_metadata.R"))
hierarchy <- get_location_metadata(location_set_id = 111, location_set_version_id = 771, release_id = 9)
us_terr <- c("American Samoa", "Guam", "Northern Mariana Islands", "Puerto Rico", "United States Virgin Islands")
midas_locs <- fread("https://raw.githubusercontent.com/midas-network/covid19-scenario-modeling-hub/master/data-locations/locations.csv")

# scenarios <- data.table(expand.grid(c("mod_npi", "low_npi"), c("high_vacc", "low_vacc")))
# scenarios[, scen := paste0(Var1, "_", Var2)]
# scenarios <- rbind(scenarios, data.table(scen = "reference"), fill = T)
# scenarios$scenario_name <- c("highVac_modNPI", "highVac_lowNPI", "lowVac_modNPI", "lowVac_lowNPI", "reference")
# scenarios$scenario_id <- c("A-2021-05-02", "B-2021-05-02", "C-2021-05-02", "D-2021-05-02", NA)

scenarios <- data.table(expand.grid(c("mod_npi"), c("low_vacc")))
scenarios[, scen := paste0(Var1, "_", Var2)]
scenarios <- rbind(data.table(scen = "reference"), scenarios, fill = T)
scenarios$scenario_name <- c("lowVac_modNPI", "lowVac_lowNPI")
scenarios$scenario_id <- c( "C-2021-05-02", "D-2021-05-02")

measures <- c("deaths", "cases", "hospitalizations")
file_names <- c("unscaled_daily_deaths.csv", "daily_infected.csv", "hospital_admissions.csv")
names(file_names) <- measures

week_dates <- seq(start_date - 1, end_date, by = "1 week")
week_dt <- data.table(week_start = head(week_dates, -1), week_end = tail(week_dates, -1))
week_dt[, week_num := .I]
week_date_dt <- rbindlist(lapply(seq(start_date, end_date, by = "1 day"), function(d) {
  copy(week_dt)[d > week_start & d <= week_end][, date := d]
}))

idr <- fread(file.path(root_dir, output_version, "reference/output_summaries/infection_detection_ratio.csv"))
idr[, date := as.Date(date)]
setnames(idr, "mean", "idr")

# Read files
dt <- rbindlist(lapply(scenarios$scen, function(s) {
  scen_dt <- rbindlist(lapply(file_names, function(f) {
    path <- file.path(root_dir, output_version, s, "output_draws", f)
    measure_dt <- fread(path)
    measure_dt[, measure := names(file_names[which(file_names == f)])]
  }), fill = T)
  scen_dt[, scen := s]
}))
dt[, date := as.Date(date)]

merge_dt <- merge(dt, idr[, .(location_id, date, idr)], by = c("location_id","date"))
merge_dt[is.na(idr), idr := 1]
merge_dt[measure == "cases", paste0("draw_", 0:99) := 
  merge_dt[measure == "cases", paste0("draw_", 0:99)] * merge_dt[measure == "cases"]$idr]
merge_dt[, idr := NULL]

scen_name_dt <- merge(merge_dt, scenarios[, .(scen, scenario_name, scenario_id)], by = c("scen"))
week_name_dt <- merge(scen_name_dt, week_date_dt, by = "date")
week_results <- week_name_dt[, lapply(.SD, sum), 
  by = .(scenario_name, scenario_id, measure, location_id, week_end, week_num),
  .SDcols = paste0("draw_", 0:99)]

week_results[measure == "deaths", target := paste0(week_num, " wk ahead inc death")]
week_results[measure == "cases", target := paste0(week_num, " wk ahead inc case")]
week_results[measure == "hospitalizations", target := paste0(week_num, " wk ahead inc hosp")]
setnames(week_results, "week_end", "target_end_date")
week_results[, week_num := NULL]

draw_mat <- t(as.matrix(week_results[, paste0("draw_", 0:99)]))
quant_list <- mclapply(quants, function(q) {
  print(q)
  apply(draw_mat, 2, quantile, q)
}, mc.cores = 10)
quant_cols <- do.call(cbind, quant_list)
colnames(quant_cols) <- format(round(quants, 3), digits = 3)
mean_col <- apply(draw_mat, 2, mean)
summ_stats_cols <- cbind(quant_cols, mean_col)

summ_dt <- as.data.table(cbind(
  week_results[,setdiff(names(week_results), paste0("draw_", 0:99)), with = F],
  summ_stats_cols
))

melt_week <- melt(summ_dt, id.vars = setdiff(names(week_results), paste0("draw_", 0:99)), variable.name = "quantile")
melt_week[quantile == "mean_col", c("quantile", "type") := .(NA, "point")]
melt_week[is.na(type), type := "quantile"]
# Plot
loc_list <- intersect(
  hierarchy[location_id == 102 | parent_id == 102 | location_name %in% us_terr]$location_id,
  unique(melt_week$location_id)
)
pdf("/homes/aucarter/Downloads/cdc_scenario.pdf", width = 11, height = 8)
for(loc_id in loc_list) {
  print(loc_id)
  plot_dt <- melt_week[location_id == loc_id & type == "point"]
  gg <- ggplot(plot_dt, aes(x = target_end_date, y = value, color = scenario_name)) + geom_line() +
    facet_wrap(~measure, scales = "free_y", nrow = 3) +
    geom_vline(xintercept = as.Date("2021-10-30"), alpha = 0.5) + 
    theme_bw() +
    theme(legend.position="bottom", legend.title = element_blank()) + 
    ggtitle(hierarchy[location_id == loc_id]$location_name) +
    xlab("Date") + ylab("Weekly counts") + scale_y_continuous(limits=c(0, NA))
  print(gg)
}
dev.off()

name_dt <- merge(melt_week[location_id %in% loc_list], hierarchy[, .(location_id, location_name)])
# setdiff(midas_locs$location_name, unique(name_dt$location_name))
# setdiff(unique(name_dt$location_name), midas_locs$location_name)
name_dt[location_name == "United States of America", location_name := "US"]
name_dt[location_name == "United States Virgin Islands", location_name := "Virgin Islands"]
midas_name_dt <- merge(name_dt, midas_locs[, .(location, location_name)], by = "location_name", all.x = T)
midas_name_dt[, c("location_name", "location_id") := NULL]

midas_name_dt[, model_projection_date := start_date]
out_dt <- midas_name_dt[measure %in% c("cases", "deaths") & scenario_name != "reference", .(model_projection_date, target, target_end_date, scenario_name, scenario_id, location, type, quantile, value)]

out_dir <- "/homes/aucarter/projects/covid/covid19-scenario-modeling-hub/data-processed/IHME-IHME COVID model (deaths unscaled)"
out_file <- paste0(start_date, "-IHME-IHME COVID model (deaths unscaled).csv")
out_path <- file.path(out_dir, out_file)
out_dt[location %in% 1:9, location := paste0("0", location)]
write.csv(out_dt, out_path, row.names = F)
