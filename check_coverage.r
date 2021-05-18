library(data.table); library(ggplot2); library(parallel)
root_dir <- "/ihme/covid-19-2/seir-covariates/latest/vaccine_coverage"


high_dt <- fread(file.path(root_dir,"vaccinations_fast_balanced_info.csv"))
high_dt[, scenario := "fast_balanced"]
low_dt <- fread(file.path(root_dir,"vaccinations_reference_info.csv"))
low_dt[, scenario := "reference"]

dt <- rbind(high_dt, low_dt)
dt[, date := as.Date(date)]
dt[, value := NULL]

melt_dt <- melt(dt, id.vars = c("location_id", "date", "scenario"), value.var = "value")

plot_dt <- melt_dt[location_id == 102]
plot_dt <- plot_dt[, .(value = sum(value)), by = .(date, scenario)]
plot_dt[, cum_value := cumsum(value), by = .(scenario)]
gg <- ggplot(plot_dt, aes(x = date, y = cum_value, color = scenario)) + geom_line()
  # facet_wrap(~variable, scales = "free_y", nrow = 5) + theme(legend.position = "bottom")
print(gg)
dev.off()
