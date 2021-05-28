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
us_dt <- vaccine_dt[location_id %in% 523:573, c("location_id", "date", "over65_population", "under65_population", grep(".*r_vaccinated", names(vaccine_dt), value = T)), with = F]
# TODO: This needs to be augmented to scale up to 95% in each bin, retaining state proportions, to bring the national coverage to 83%
dt <- as.data.table(rbind(
  us_dt[, .(location_id, date, under65_population, lr_vaccinated)] %>%
    setnames(., c("under65_population", "lr_vaccinated"),
      c("population", "protected")) %>%
      mutate(group = "lr"),
  us_dt[, .(location_id, date, over65_population, hr_vaccinated)] %>%
    setnames(., c("over65_population", "hr_vaccinated"),
      c("population", "protected")) %>%
      mutate(group = "hr")
))
dt[, cum_protected := cumsum(protected), by = .(location_id, group)]
dt[, coverage := cum_protected / population]

logit <- function(x) log(x / (1 - x))
ilogit <- function(x) exp(x) / (1 + exp(x))

obj_fun <- function(k, props, weights, target, cap) {
  shifted_props <- ilogit(logit(props) + k) * cap 
  w_mean <- weighted.mean(shifted_props, weights)
  return((target - w_mean)**2)
}

logit_shift <- function(props, k, cap) {
  ilogit(logit(props) + k) * cap
}

logit_rake <- function(props, weights, target, cap) {
  sol <- optimize(
    obj_fun, 
    interval = c(-1, 1) * 1e2,
    props = props, weights = weights, target = target, cap = cap
  )
  k <- sol$minimum
  opt_shifted_props <- logit_shift(props, k, cap)
  return(list(opt_shifted_props, k))
}

date_dt <- dt[date == "2021-10-31"]
target = 0.83
cap = 0.95
raked_props <- logit_rake(date_dt$coverage, date_dt$population, target, cap)

weighted.mean(raked_props[[1]], date_dt$population)

dt[, raked_coverage := logit_shift(coverage, raked_props[[2]], cap)]
dt[, date := as.Date(date)]
dt[, c("protected", "cum_protected") := NULL]
melt_dt <- melt(dt, id.vars = c("location_id", "date", "group", "population"))

date_dt <- dt[date == "2021-10-31"]
weighted.mean(date_dt$raked_coverage, date_dt$population)

library(ggplot2)
plot_dt <- melt_dt[location_id == 555]
pdf("/homes/aucarter/Downloads/raked_coverage.pdf", height = 20, width = 9)
gg <- ggplot(melt_dt, aes(x = date, y = value, linetype = group, color = variable)) + 
  geom_line() + facet_wrap(~location_id, nrow = 10) + 
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0.95, col = "red") + ylim(c(0, 1))
print(gg)
dev.off()
