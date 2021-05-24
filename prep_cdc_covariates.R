library(data.table)

## Scenario parameters
scen_eff_reduction <- 0.2
scen_eff_label <- "cdc_low"
scen_eff_reduction_start <- as.Date("2021-04-30")
scen_eff_reduction_end <- as.Date("2021-11-01")

reduction_dt <- as.data.table(
  approx(
    x = c(scen_eff_reduction_start, scen_eff_reduction_end), 
    y = c(1, scen_eff_reduction),
    xout = seq(scen_eff_reduction_start,scen_eff_reduction_end, by = "1 day")
  )
)
names(reduction_dt) <- c("date", "reduction")


## Prep mobility
mob_path <- "/ihme/covid-19-2/mobility-covariate/best/mobility_reference.csv"
mob_dt <- fread(mob_path)
mob_dt[, date := as.Date(date)]

scen_mob <- merge(mob_dt, reduction_dt, by = "date", all.x = T)
scen_mob[date > scen_eff_reduction_end, reduction := scen_eff_reduction]
scen_mob[date > scen_eff_reduction_end, reduction := scen_eff_reduction]
eff_vars <- c("mobility_forecast")
scen_mob[!is.na(reduction), (eff_vars) := 
  scen_mob[!is.na(reduction), eff_vars, with = F] * reduction]
scen_mob[, reduction := NULL]
mob_out_path <- paste0("/ihme/covid-19-2/mobility-covariate/best/mobility_", scen_eff_label,".csv")
write.csv(scen_mob, mob_out_path, row.names = F)

## Prep mask use
mask_path <- "/ihme/covid-19-2/mask-use-outputs/best/mask_use.csv"
mask_dt <- fread(mask_path)
mask_dt[, date := as.Date(date)]

scen_mask <- merge(mask_dt, reduction_dt, by = "date", all.x = T)
scen_mask[date > scen_eff_reduction_end, reduction := scen_eff_reduction]
scen_mask[date > scen_eff_reduction_end, reduction := scen_eff_reduction]
eff_vars <- c("mask_use")
scen_mask[!is.na(reduction), (eff_vars) := 
  scen_mask[!is.na(reduction), eff_vars, with = F] * reduction]
scen_mask[, reduction := NULL]
mask_out_path <- paste0("/ihme/covid-19-2/mask-use-outputs/best/mask_use_", scen_eff_label,".csv")
write.csv(scen_mask, mask_out_path, row.names = F)

## Prep vaccine coverage
vaccine_coverage_path <- "/ihme/covid-19-2/vaccine-coverage/best/slow_scenario_vaccine_coverage.csv"
vaccine_dt <- fread(vaccine_coverage_path)
names(vaccine_dt)
us_dt <- vaccine_dt[location_id %in% 523:573, c("location_id", "date", grep("_effective_protected", names(vaccine_dt), value = T)), with = F]
# TODO: This needs to be augmented to scale up to 95% in each bin, retaining state proportions, to bring the national coverage to 83%