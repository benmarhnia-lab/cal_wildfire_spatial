##################
# This document includes codes written by Chen Chen for analyses in the manuscript
# "Spatial Heterogeneity of the Respiratory Health Impacts of Wildfire Smoke PM2.5 in California"
# written by Chen Chen on 08/21/23
##################

indir1 <- "" ## directory for exposure data and community characteristic data
indir2 <- "" ## directory for outcome data
outdir1 <- "~/cal_wildfire_spatial" ## directory to store results
outdir2 <- file.path(outdir1, "figures", "spatial") ## directory to store zcta-specific results
if (!dir.exists(outdir2)) dir.create(outdir2)

dataset <- "wf15_binary_1773zcta"

library(data.table)
library(sf)
library(survival)
library(msm)
library(splines)

## run case-crossover analysis at state-level
##################
if (!file.exists(file.path(outdir1, "results", "state_models"))) dir.create(file.path(outdir1, "results", "state_models"))

## read in exposure dataset
exposure <- readRDS(file.path(outdir1, "data", paste0(dataset, "_0619.rds")))

## read in the list of zctas to include based on my standard
all <- fread(file.path(outdir1, "results", paste0("zcta_list_pop1000_someexposure_", dataset, ".csv")))[, .(zcta, pop)] ## main analysis
exposure <- exposure[exposure$zcta %in% all$zcta, ] ## only include zctas within the list for analysis

dataset <- paste0(dataset, "_", length(unique(exposure$zcta)))

## clean health data to create dataset for case-crossover
ha <- fread(file.path(indir2, "combo_99_19_patzip_cvd_resp_EM.csv"))
names(ha)[1:2] <- c("zcta", "case_date")
## Focus on overall population
ha <- ha[ha$case_date > as.Date("2005-12-31") & ha$case_date < as.Date("2020-01-01"), .(zcta, case_date, respiratory)]
length(unique(ha$zcta)) ## 1778 zipcodes
ha <- ha[ha$zcta %in% exposure$zcta, ] ## keep of data from zipcodes with exposure data
length(unique(ha$zcta)) ## 1396 zipcodes

## create variables for merging
ha$wday <- wday(ha$case_date)
ha$month <- month(ha$case_date)
ha$year <- year(ha$case_date)

## create empty variable 
out <- numeric()

lags <- c("same_day", "lag1")
outcomes <- c("respiratory")
for (lag in lags) {
  if (lag!="same_day") { ## create lags in exposure
    nn <- as.numeric(gsub("lag", "", lag))
    new <- copy(exposure)
    new <- new[order(new$date), ]
    new[, wf:=shift(wf, n = nn, fill = NA, type = "lag"), by=.(zcta)]
  } else {
    new <- copy(exposure)
    new <- new[order(new$date), ]
  }
  new$wday <- wday(new$date)
  for (outcome in outcomes) {

    resp <- ha[respiratory>0 & !is.na(respiratory), ]
    resp$id <- 1:nrow(resp)
    dt <- merge(resp, new, by = c("zcta", "wday", "month", "year"), all.x = TRUE, allow.cartesian=TRUE)
    
    dt$case <- ifelse(dt$case_date==dt$date, 1, 0)
    loc <- grep(outcome, names(dt))
    names(dt)[loc] <- "outcome"
    
    m.wf <- clogit(case ~ wf + strata(id), data=dt, weights=outcome, method="approximate")
    saveRDS(m.wf, file.path(outdir1, "results", "state_models", paste0(dataset, "_", lag, "_", outcome, "_wf.rds")))
    
    m.wf.ln <- clogit(case ~ wf + tmpt + strata(id), data=dt, weights=outcome, method="approximate")
    saveRDS(m.wf.ln, file.path(outdir1, "results", "state_models", paste0(dataset, "_", lag, "_", outcome, "_wf_tmpt.rds")))
    
    m.wf.ns <- clogit(case ~ wf + ns(tmpt, 6) + strata(id), data=dt, weights=outcome, method="approximate")
    saveRDS(m.wf.ns, file.path(outdir1, "results", "state_models", paste0(dataset, "_", lag, "_", outcome, "_wf_ns6tmpt.rds")))
    
    ## create results for plotting
    temp <- data.frame(rbind(summary(m.wf)$conf.int[1, c(1, 3, 4)],
                             summary(m.wf.ln)$conf.int[1, c(1, 3, 4)],
                             summary(m.wf.ns)$conf.int[1, c(1, 3, 4)]))
    row.names(temp) <- c("wildfire", "wildfire_tmpt", "wildfire_ns6tmpt")
    names(temp) <- c("est", "ll", "ul")
    temp$se <- c(sqrt(m.wf$var[1,1]), sqrt(m.wf.ln$var[1,1]), sqrt(m.wf.ns$var[1,1]))
    temp$estimate <- row.names(temp)
    out <- rbind(out, cbind(outcome = outcome, lag = lag, temp))
    
    m.wf <- m.wf.ln <- m.wf.ns <- dt <- temp <- NULL
  }
}

write.csv(out, file.path(outdir1, "results", paste0(dataset, "_state_model_summary.csv")), row.names = FALSE)

##################

## run case-crossover analysis at county level and air basin level
##################
if (!file.exists(file.path(outdir1, "results", "region_models"))) dir.create(file.path(outdir1, "results", "region_models"))

## read in the list of zctas to include based on my standard
all <- fread(file.path(outdir1, "results", paste0("zcta_list_pop1000_someexposure_", dataset, ".csv")))[, .(zcta, pop)] ## main analysis

## Create connection between ZCTA and county based on ZCTA5 county relationship file downloaded from U.S census 2010
county <- fread(file.path(indir1, "contours","ZIP_County_CA.csv"))
ref <- fread(file.path(indir1, "contours", "www2.census.gov_geo_docs_reference_codes_files_st06_ca_cou.txt"), 
             colClasses=list(character=2:3))
ref$fips <- paste0(ref$V2, ref$V3)
ref$county <- gsub(" County", "", ref$V4)
county <- merge(county, ref[, c("county", "fips")], by="county", all.x=TRUE)
connect <- county[, c("zip", "county", "fips")]
names(connect) <- c("zcta", "region", "code")

## create connection between ZCTA and air basin
air.basin <- as.data.frame(read_sf(file.path(indir1, "contours", "CA_ZIP_PWC_AirBasin",
                                             "CA_ZIP_PWC_AirBasin.shp"), stringsAsFactors = FALSE, as_tibble = FALSE))[, c("Join_Count", "ZCTA5", "NAME", "CAABA_ID")]
names(air.basin)[2:4] <- c("zcta", "region", "code")
print(air.basin [air.basin$zcta %in% c(93920, 94590, 93924, 94592),]) ## first two zipcodes have no air basin data so I used info from nearby zipcodes based on google map
air.basin[air.basin$zcta==93920, c("region", "code")] <- air.basin[air.basin$zcta==93924, c("region", "code")]
air.basin[air.basin$zcta==94592, c("region", "code")] <- air.basin[air.basin$zcta==94590, c("region", "code")]
connect <- rbind(connect, air.basin[, c("zcta", "region", "code")])

connect <- connect[connect$zcta %in% all$zcta, ]
regions <- unique(connect$region) ## note "alpine"/06003 county not included because both zipcodes in this county were excluded

## read in exposure dataset
exposure <- readRDS(file.path(outdir1, "data", paste0(dataset, "_0619.rds")))
exposure <- exposure[exposure$zcta %in% all$zcta, ] ## only include zctas within the list for analysis

dataset <- paste0(dataset, "_", length(unique(exposure$zcta)))

## clean health data to create dataset for case-crossover
ha <- fread(file.path(indir2, "combo_99_19_patzip_cvd_resp_EM.csv"))
names(ha)[1:2] <- c("zcta", "case_date")
## Focus on overall population
ha <- ha[ha$case_date > as.Date("2005-12-31") & ha$case_date < as.Date("2020-01-01"), .(zcta, case_date, respiratory)]
length(unique(ha$zcta)) ## 1778 zipcodes
ha <- ha[ha$zcta %in% exposure$zcta, ] ## keep of data from zipcodes with exposure data
length(unique(ha$zcta)) ## 1772 zipcodes

## create variables for merging
ha$wday <- wday(ha$case_date)
ha$month <- month(ha$case_date)
ha$year <- year(ha$case_date)

## create empty variable 
out <- numeric()

lags <- c("same_day", "lag1")
outcomes <- c("respiratory")
for (lag in lags) {
  if (lag!="same_day") { ## create lags in exposure
    nn <- as.numeric(gsub("lag", "", lag))
    new <- copy(exposure)
    new <- new[order(new$date), ]
    new[, wf:=shift(wf, n = nn, fill = NA, type = "lag"), by=.(zcta)]
  } else {
    new <- copy(exposure)
    new <- new[order(new$date), ]
  }
  new$wday <- wday(new$date)
  for (outcome in outcomes) {
    
    resp <- ha[respiratory>0 & !is.na(respiratory), ]
    resp$id <- 1:nrow(resp)
    dt <- merge(resp, new, by = c("zcta", "wday", "month", "year"), all.x = TRUE, allow.cartesian=TRUE)
    
    dt$case <- ifelse(dt$case_date==dt$date, 1, 0)
    loc <- grep(outcome, names(dt))
    names(dt)[loc] <- "outcome"
    
    for (region_ in regions) {
      foo <- dt[dt$zcta %in% connect$zcta[connect$region==region_], ]
      
      m.wf <- clogit(case ~ wf + strata(id), data=foo, weights=outcome, method="approximate")
      saveRDS(m.wf, file.path(outdir1, "results", "region_models", paste0(dataset, "_", lag, "_", outcome, "_", region_, "_wf.rds")))
      
      ## create results for plotting
      temp <- data.frame(rbind(summary(m.wf)$conf.int[1, c(1, 3, 4)]))
      row.names(temp) <- c("wildfire")
      names(temp) <- c("est", "ll", "ul")
      temp$se <- c(sqrt(m.wf$var[1,1]))
      temp$estimate <- row.names(temp)
      out <- rbind(out, cbind(outcome = outcome, lag = lag,
                              temp, 
                              region = region_, code = unique(connect$code[connect$region==region_])))
      
      m.wf <- m.wf.ln <- m.wf.ns <- foo <- temp <- NULL
    }
  }
}

write.csv(out, file.path(outdir1, "results", paste0(dataset, "_region_model_summary.csv")), row.names = FALSE)

##################

## attributable number of hospital admission for wildfire smoke using 
## AN = AF * sum of HA in exposed day = (1-1/RR) * sum of HA in wildfire days of all 1773 zctas
##################
out <- fread(file.path(outdir1, "results", paste0(dataset, "_1396_region_model_summary.csv")), colClasses = c("code" = "character"))
out.st <- fread(file.path(outdir1, "results", paste0(dataset, "_1396_state_model_summary.csv")))

## read in the list of zctas to include based on my standard
all <- fread(file.path(outdir1, "results", paste0("zcta_list_", dataset, ".csv")))[, .(zcta, pop)] ## starting point of 1773 zctas

## Create connection between ZCTA and county based on ZCTA5 county relationship file downloaded from U.S census 2010
county <- fread(file.path(indir1, "contours", "ZIP_County_CA.csv"))
ref <- fread(file.path(indir1, "contours", "www2.census.gov_geo_docs_reference_codes_files_st06_ca_cou.txt"), 
             colClasses=list(character=2:3))
ref$fips <- paste0(ref$V2, ref$V3)
ref$county <- gsub(" County", "", ref$V4)
county <- merge(county, ref[, c("county", "fips")], by="county", all.x=TRUE)
connect <- county[, c("zip", "county", "fips")]
names(connect) <- c("zcta", "region", "code")

## create connection between ZCTA and air basin
air.basin <- as.data.frame(read_sf(file.path(indir1, "contours", "CA_ZIP_PWC_AirBasin",
                                             "CA_ZIP_PWC_AirBasin.shp"), stringsAsFactors = FALSE, as_tibble = FALSE))[, c("Join_Count", "ZCTA5", "NAME", "CAABA_ID")]
names(air.basin)[2:4] <- c("zcta", "region", "code")
print(air.basin [air.basin$zcta %in% c(93920, 94590, 93924, 94592),]) ## first two zipcodes have no air basin data so I used info from nearby zipcodes based on google map
air.basin[air.basin$zcta==93920, c("region", "code")] <- air.basin[air.basin$zcta==93924, c("region", "code")]
air.basin[air.basin$zcta==94592, c("region", "code")] <- air.basin[air.basin$zcta==94590, c("region", "code")]
connect <- rbind(connect, air.basin[, c("zcta", "region", "code")])

connect <- connect[connect$zcta %in% all$zcta, ]
regions <- unique(connect$region) ## note "alpine"/06003 county not included because both zipcodes in this county were excluded

## read in exposure dataset
exposure <- readRDS(file.path(outdir1, "data", paste0(dataset, "_0619.rds")))
exposure <- exposure[exposure$zcta %in% all$zcta, .(zcta, date, wf)] ## only include zctas within the list for analysis

## clean health data to create dataset for case-crossover
ha <- fread(file.path(indir2, "combo_99_19_patzip_cvd_resp_EM.csv"))
names(ha)[1:2] <- c("zcta", "date")
## Focus on overall population
ha <- ha[ha$date > as.Date("2005-12-31") & ha$date < as.Date("2020-01-01"), .(zcta, date, respiratory)]
length(unique(ha$zcta)) ## 1778 zipcodes
ha <- ha[ha$zcta %in% exposure$zcta, ]

lags <- c("same_day", "lag1")
outcomes <- c("respiratory")
out$exposed.case <- 0
out.st$exposed.case <- 0
for (lag in lags) {
  if (lag!="same_day") { ## create lags in exposure
    nn <- as.numeric(gsub("lag", "", lag))
    new <- copy(exposure)
    new <- new[order(new$date), ]
    new[, wf:=shift(wf, n = nn, fill = NA, type = "lag"), by=.(zcta)]
  } else {
    new <- copy(exposure)
    new <- new[order(new$date), ]
  }
  dt <- merge(ha, new, by = c("zcta", "date"), all.x = TRUE, allow.cartesian=TRUE)
  
  for (outcome in outcomes) {
    loc.st <- which(out.st$outcome==outcome & out.st$lag==lag)
    out.st[loc.st, "exposed.case"] <- dt[dt$wf == 1, sum(eval(as.name(outcome)))]
    
    for (region_ in regions) {
      foo <- dt[dt$zcta %in% connect$zcta[connect$region==region_] & dt$wf == 1, sum(eval(as.name(outcome)))]
      loc <- which(out$outcome==outcome & out$lag==lag & out$region==region_)
      out[loc, "exposed.case"] <- foo
    }
  }
}

out$paf <- 1-1/out$est
out$paf.ll <- 1-1/out$ll
out$paf.ul <- 1-1/out$ul

out$pan <- out$paf * out$exposed.case
out$pan.ll <- out$paf.ll * out$exposed.case
out$pan.ul <- out$paf.ul * out$exposed.case

out.st$paf <- 1-1/out.st$est
out.st$paf.ll <- 1-1/out.st$ll
out.st$paf.ul <- 1-1/out.st$ul

out.st$pan <- out.st$paf * out.st$exposed.case
out.st$pan.ll <- out.st$paf.ll * out.st$exposed.case
out.st$pan.ul <- out.st$paf.ul * out.st$exposed.case

write.csv(out, file.path(outdir1, "results", paste0(dataset, "_1396rr_1773pan_region_model_summary.csv")), row.names = FALSE)
write.csv(out.st, file.path(outdir1, "results", paste0(dataset, "_1396rr_1773pan_state_model_summary.csv")), row.names = FALSE)
##################

## run within-community matched design (Schwarz et al. 2021) after incorporating
## control exclusion criteria from Bobb et al. 2014 and Liu et al. 2017
## Schwarz, L., Hansen, K., Alari, A., Ilango, S. D., Bernal, N., Basu, R., et al. (2021). Spatial variation in the joint effect of extreme heat events and ozone on respiratory hospitalizations in California. Proceedings of the National Academy of Sciences, 118(22). https://doi.org/10.1073/pnas.2023078118
## Bobb, J. F., Obermeyer, Z., Wang, Y., & Dominici, F. (2014). Cause-Specific Risk of Hospital Admission Related to Extreme Heat in Older Adults. JAMA, 312(24), 2659–2667. https://doi.org/10.1001/jama.2014.15715
## Liu, J. C., Wilson, A., Mickley, L. J., Dominici, F., Ebisu, K., Wang, Y., et al. (2017). Wildfire-specific Fine Particulate Matter and Risk of Hospital Admissions in Urban and Rural Counties: Epidemiology, 28(1), 77–85. https://doi.org/10.1097/EDE.0000000000000556
##################
## identify potential control days for each exposed day baed on schwartz et al. and Bobb et al. method: controls identified as 
## 1) within the window of buffer calendar days (nbuffer) before or after the exposed day and 
## 2) separated from any other exposed day for more than 3 days
month.control <- function(exposed, control, nbuffer) { 
  if (length(exposed) > 0) {
    out <- lapply(exposed, function(e_day) {
      e_buffer <- e_day + (-nbuffer:nbuffer) ## potential days for control (buffer day window from exposure)
      c_pool <- control[control %in% e_buffer]  ## potential control days after excluding days close to other exposure
      return(c(e_day, c_pool))
    })
    names(out) <- exposed
  } else {
    out <- numeric()
  }
  return(out)
}

## weight the outcome of all potential controls using the distance in day and 
## directly calculate the incidence ratio for each match exposed and control group, then calculate the average
## updated error message so that zctas with exposed day but no eligible control day would be marker "no available control day" on 1/11/23 -- no trial run yet
month.wt.analysis <- function(exposure, event, c.list, outcome.dt) { 
  if (length(c.list) > 0) {
    out <- numeric()
    for (i in 1:length(c.list)) {
      if (length(c.list[[i]]) > 1) {
        e_day <- c.list[[i]][1]
        c_pool <- data.table(date=c.list[[i]], exposed = c(1, rep(0, length(c.list[[i]])-1)))
        c_pool$wt <- 1/round(abs(as.numeric(difftime(c_pool$date, e_day, units = "days"))), digits = 0)
        c_pool$wt[1] <- 0 ## assign 0 weight to the day with exposure
        c_pool <- merge(c_pool, outcome.dt, by="date", all.x=TRUE)
        
        rd <- c_pool[exposed==1, eval(as.name(event))] - c_pool[, sum(eval(as.name(event))*wt)/sum(wt)]
        rd <- ifelse(is.infinite(rd), NA, rd)
        out <- c(out, rd)
        
      }
    }
    
    if (length(out)==0) {
      temp <- simpleError(paste("No available control day"))
    } else {
      temp <- c(mean(out, na.rm=TRUE), sum(!is.na(out)))
      if (is.na(temp[1])) {
        temp <- simpleError(paste("No non-infinite RR"))
        temp$call <- NULL
      }
    }
    
    
  } else {
    temp <- simpleError(paste("No available exposed day"))
    temp$call <- NULL
  }
  return(temp)
}

dt <- readRDS(file.path(outdir1, "data", paste0(dataset, "_0619.rds")))
dt$yday <- yday(dt$date)

## Focus on overall population
ha <- fread(file.path(indir2, "combo_99_19_patzip_cvd_resp_EM.csv"))
names(ha)[1:2] <- c("zcta", "date")
ha <- ha[ha$date > as.Date("2005-12-31") & ha$date < as.Date("2020-01-01"), .(zcta, date, respiratory)]

bobb.control.n <- 4

events <- c("respiratory")
nms <- c(paste(rep(c("month_wt"), each=2), rep(c("rd", "ngrp"), times=1), sep="_"))
nms <- paste(rep(nms, times=1), rep(c("wf1"), each=length(nms)), sep="_")
bar <- data.frame(zcta=rep(unique(dt$zcta), each = length(events)), event = events, 
                  wf1= NA, 
                  bobb_control_pool = NA, 
                  setNames(replicate(length(nms), NA, simplify = F), nms)
)
fail <- numeric()
set.seed(824)
for (i in 1:nrow(bar)) {
  event_ <- bar$event[i]
  ha_ <- ha[ha$zcta==bar$zcta[i], ]
  foo <- dt[dt$zcta==bar$zcta[i], ]
  foo <- merge(foo, ha_[, .(date, respiratory)], by="date", all.x = TRUE)
  setnafill(foo, fill=0, cols = c("respiratory")) ## fill in zeros for days without resp

  ## potential control identification method from Bobb et al--remove those close to exposure
  buffer_days <- foo$date[foo$wf==1]
  buffer_days <- unique(c(buffer_days-3, buffer_days-2, buffer_days-1, buffer_days, buffer_days + 1, buffer_days + 2, buffer_days + 3))
  day00_bobb <- foo$date[!(foo$date %in% buffer_days)]
  bar$bobb_control_pool[i] <- length(day00_bobb) ## total number of potential control days excluding those close to exposure
  
  for (exposure_ in c("wf1")) {
    ## identify exposed days using Bobb method
    days <- foo$date[foo$wf==as.numeric(substring(exposure_, 3, 3))]
    bar[i, exposure_] <- length(days)
    
    ## control identification method based on Schwarz et al. 2021: same year within 30 days before and after the exposure and exclude these within three days of an event
    baz2 <- month.control(days, day00_bobb, 30)
    
    ## weighted analysis method adopted from Schwarz et al. 2021
    temp4 <- month.wt.analysis(exposure_, event_, baz2, foo)
    if (!inherits(temp4, what = "condition")) {
      bar[i, paste0(c("month_wt_rd_", "month_wt_ngrp_"), exposure_)] <- temp4
    } else {
      fail <- rbind(fail, data.frame(zcta=bar$zcta[i], event = event_, exposure = exposure_, ngrps = length(baz2), method = "month_wt", error=as.character(temp4)))
      bar[i, paste0("month_wt_ngrp_", exposure_)] <- length(baz2)
    }
    
  }
  
  if (i%%500==0) cat("\n", "first two methods finished zcta #", i, "\n")
}

write.csv(bar, file.path(outdir1, "results", paste0(dataset, "_specific.csv")), row.names = FALSE)
write.csv(fail, file.path(outdir1, "results", paste0(dataset, "_specific_fail.csv")), row.names = FALSE)
##################


## Apply Bayesian spatial analysis to incorporate spatial heterogeneity into the estimates for risk difference
## Schwarz, L., Hansen, K., Alari, A., Ilango, S. D., Bernal, N., Basu, R., et al. (2021). Spatial variation in the joint effect of extreme heat events and ozone on respiratory hospitalizations in California. Proceedings of the National Academy of Sciences, 118(22). https://doi.org/10.1073/pnas.2023078118
library(sp)
library(gstat) 
library("spBayes")
library(coda)
library(MBA)

if (!dir.exists(file.path(outdir2, dataset))) dir.create(file.path(outdir2, dataset))
##################
events <- c("respiratory")

methods <- c("month_wt")

## read in data for analysis
bar <- fread(file.path(outdir1, "results", paste0(dataset, "_specific.csv")))
fail <- fread(file.path(outdir1, "results", paste0(dataset, "_specific_fail.csv")))
coords <- read.csv(file.path(outdir1, "results", paste0("zcta_list_", dataset, ".csv")), as.is = TRUE) ## same as above but with population info
bar <- merge(bar, coords, by = "zcta", all.x=TRUE)

for (event_ in events) {
  outdir3 <- file.path(outdir2, dataset, event_)
  if (!dir.exists(outdir3)) dir.create(outdir3)
  sink(file.path(outdir3, "summary of running spatial bayesian_flatprior.txt"))
  for (m in methods) {
    set.seed(824)
    cat("\n\n", event_, m, "\n")
    fail_m <- fail[fail$method==m & fail$event==event_, ]
    fail.zcta <- unique(fail_m$zcta[fail_m$ngrps==0])
    baz <- bar[!(zcta %in% fail.zcta) & event==event_, ]
    cat("# of zctas removed due to no exposed day is", length(fail.zcta), "\n")
    nn <- nrow(baz)
    baz <- baz[baz$pop>1000 &!is.na(baz$pop), ]
    cat("# of zctas removed due to no pop or pop<=1000 is", nn-nrow(baz), "\n")
    cat("# of zctas left after 1k pop and exposure day criteria is", nrow(baz), "\n")
    ex.fail <- fail_m[fail_m$ngrps>0]
    ex.fail <- ex.fail[!ex.fail$zcta %in% fail.zcta, ]
    ex.fail <- ex.fail[ex.fail$zcta %in% baz$zcta, ]
    cat("# of zctas removed due to failures other than no exposed day or pop is", length(unique(ex.fail$zcta)), "\n")
    cat("zctas removed due to failures other than no exposed day or pop", "\n")
    print(merge(ex.fail, coords[,c("zcta", "pop")], by="zcta", all.x=TRUE))
    baz <- baz[!baz$zcta %in% ex.fail$zcta, ]
    
    ## calculate rate difference from count difference, per 100000 person
    baz[, paste0(m, "_srd"):= eval(as.name(paste0(m, "_rd_wf1")))/pop * 100000]
    
    cat("# of zctas with estimates for", m, "is", nrow(baz), "\n")
    
    bayesDF <- baz[, c("zcta", paste0(m, "_srd"), "FinalLon", "FinalLat"), with=FALSE]
    names(bayesDF)[2] <- "est"
    
    # US census bureau datum
    spdf <- SpatialPointsDataFrame(coords = bayesDF[,.(FinalLon, FinalLat)],
                                   data = bayesDF)
    
    v1 <- variogram(est ~ 1, data = spdf)
    png(file.path(outdir3, paste0(m, nrow(bayesDF), "zcta_srd_variogram.png")))
    print(plot(v1, main=paste0(dataset, "_", m, nrow(bayesDF), "zcta"), cex=1.5))
    dev.off()
    
    ## conduct Bayesian model with flat priors
    # tau sq is nugget
    # sigma sq sill
    # phi range
    n.samples = 10000
    
    bef.sp <- spLM(est ~ 1, data = bayesDF, coords = as.matrix(bayesDF[,.(FinalLon, FinalLat)]),
                   starting = list("phi" = 2, "sigma.sq" = 16, "tau.sq" = 8),
                   tuning = list("phi" = 0.002, "sigma.sq" = 0.016, "tau.sq" = 0.008),
                   priors = list("phi.Unif" = c(0.001, 6), "sigma.sq.IG" = c(0.001, 0.001)
                                 , "tau.sq.IG" = c(0.001, 0.001)
                   ), ## flat prior for IG
                   # priors = list("phi.Unif" = c(0.001, 6), "sigma.sq.IG" = c(2, 1/16)
                   #               , "tau.sq.IG" = c(2, 1/8)
                   # ), ## actually used--(2, 1/mean) priors for IG
                   cov.model = "spherical", n.samples = n.samples, verbose = TRUE, n.report=2000)
    
    m <- paste0("flatprior_", m) ## comment out if not flatprior, important to differentiate file name
    
    cat("summary of thetas before burn-in", "\n")
    print(round(summary(mcmc(bef.sp$p.theta.samples))$quantiles, 3))
    png(file.path(outdir3, paste0(m, nrow(bayesDF), "zcta_srd_mcmctrace.png")))
    plot(bef.sp$p.theta.samples)
    dev.off()
    
    ## exclude 75% samples as burn-in
    burn.in <- floor(0.75*n.samples)
    bef.sp <- spRecover(bef.sp, start = burn.in, n.report=1000)
    
    cat("summary of thetas after burn-in", "\n")
    print(round(summary(mcmc(bef.sp$p.theta.recover.samples))$quantiles, 3))
    png(file.path(outdir3, paste0(m, nrow(bayesDF), "zcta_srd_mcmctrace_afterburnin.png")))
    plot(bef.sp$p.theta.recover.samples)
    dev.off()
    
    beta.samples <- bef.sp$p.beta.recover.samples
    w.samples <- bef.sp$p.w.recover.samples
    
    bayesDF$overall_mu <- mean(beta.samples)
    bayesDF$overall_sd <- sd(beta.samples)
    cat("Statewise mean is", mean(beta.samples), ", sd is", sd(beta.samples), "\n")
    bayesDF$w_hat_mu <- apply(w.samples, 1, mean)
    bayesDF$w_hat_sd <- apply(w.samples, 1, sd)
    bayesDF$SNR <- bayesDF$w_hat_mu/bayesDF$w_hat_sd
    bayesDF$truncSNR <- ifelse(bayesDF$SNR < 2 & bayesDF$SNR > -2, NA, bayesDF$SNR )
    print(summary(bayesDF$est))
    print(summary(bayesDF$w_hat_mu))
    print(summary(bayesDF$w_hat_sd))
    print(summary(bayesDF$SNR))
    print(summary(bayesDF$truncSNR))
    write.csv(bayesDF, file.path(outdir1, "results", paste0(dataset, "_srd_", event_, "_", m, nrow(bayesDF), "zcta.csv")))
  }
  sink()
}
##################

## Meta regression
## In summary, we transformed all indicators so that higher value represents better  
## socioeconomical status (SES) based on previous believes, except for the ethnicity groups and population density.  
library(rgdal)
library(meta)
##################
dataset <- "wf15_binary_1773zcta"
methods <- c("month_wt")
nc <- c(1396) ## 7/24/23 results
names(nc) <- methods
events <- c("respiratory")

## read in pre-pooling estimates
bar <- fread(file.path(outdir1, "results", paste0(dataset, "_specific.csv")))
coords <- read.csv(file.path(outdir1, "results", paste0("zcta_list_", dataset, ".csv")), as.is = TRUE) ## same as above but with population info
bar <- merge(bar, coords, by = "zcta", all.x=TRUE)



## read in variables for meta regression--hpi3-hpi
hpi <- fread(file.path(outdir1, "data", "zip_selected_hpi3_hpi_wt_041823.csv"))
## population density
zcta.sp <- readOGR(file.path(indir1, "contours", "tl_2010_06_zcta510",
                             "tl_2010_06_zcta510.shp"), stringsAsFactors = FALSE)
area <- zcta.sp@data[, c("ZCTA5CE10", "ALAND10")]
area$ZCTA5CE10 <- as.numeric(area$ZCTA5CE10)
hpi <- merge(hpi, area, by.x="ZIP", by.y="ZCTA5CE10", all.x=TRUE) 
hpi$ALAND10 <- as.numeric(hpi$ALAND10)
hpi$popden <- hpi$pop/hpi$ALAND10 * 10000 ## so that it's per 10k
## after transformation, all variables are true to its meaning and higher better
hpi_vbs <- c("employed", "abovepoverty", "bachelorsed", "inhighschool",
             "treecanopy", "AC", "insured", "percapitaincome", "automobile",
             "white", "black", "asian", "latino", "NativeAm","PacificIsl",
              "popden")
hpi <- hpi[!is.na(employed) & !is.na(inhighschool) & !is.na(treecanopy), c("ZIP", hpi_vbs), with=FALSE]

out <- numeric()
sink(file.path(outdir2, dataset, "summary of running meta regression srd_hpi3_hpi_wt_041823_flatprior.txt"))
for (event_ in events) {
  cat("\n\n", event_, "\n")
  outdir3 <- file.path(outdir2, dataset, event_)
  for (m in methods) {
    ## srd after pooling
    bayesDF <- fread(file.path(outdir1, "results", paste0(dataset, "_srd_", event_,  "_flatprior_", m, nc[m], "zcta.csv")))
    ntemp <- nrow(bayesDF)
    cat("\n\n", "# of zctas with", m, "based srd is", ntemp, "\n")
    bayesDF <- merge(bayesDF, hpi, by.x="zcta", by.y="ZIP", all.x=TRUE)
    bayesDF <- bayesDF[!is.na(employed), ]
    cat("# of zctas included in meta regression of", m, "based srd is", nrow(bayesDF), "\n")
    cat("# of zctas removed due to missing hpi in meta regression of", m, "based srd is", ntemp - nrow(bayesDF), "\n")
    
    nms <- paste(rep(c("meta"), each=2), rep(c("coef", "se"), time=1), sep="_")
    estimate <- data.frame(vb=hpi_vbs,
                           min=NA, q1=NA, median=NA, q3=NA, max=NA, zcta5.used=NA, 
                           event = event_,
                           setNames(replicate(length(nms), NA, simplify = FALSE), nms))
    for (i in 1:nrow(estimate)) {
      if (estimate$vb[i]=="AC") {
        cache <- copy(bayesDF)
        bayesDF <- bayesDF[!is.na(bayesDF$AC), ]
        cat("# of zctas removed due to missing AC in meta regression of", m, "based srd is", nrow(cache) - nrow(bayesDF), "\n")
      }
      estimate[i, "zcta5.used"] <- nrow(bayesDF)
      
      ## meta-regression of srd after bayesian spatial pooling
      m2 <- metagen(TE = w_hat_mu, seTE = w_hat_sd, studlab = zcta, data = bayesDF)
      g2 <- metareg(m2, as.formula(paste0(" ~ ", estimate$vb[i])))
      estimate[i, nms[3:4]] <- c(g2$beta[2], g2$se[2])
      
      if (estimate$vb[i]=="AC") {
        bayesDF <- copy(cache)
      }
    }
    estimate$iqr <- estimate$q3 - estimate$q1
    temp <- data.frame(method=m)
    temp <- cbind(rep(temp), estimate)
    out <- rbind(out, temp)
    estimate <- nms <- g <- g2 <- m2 <- temp <- NULL
  }
}
write.csv(out, file.path(outdir1, "results",  paste0(dataset, "_srd_metaregression_hpi3_hpi_wt_041823_flatprior.csv")), row.names = FALSE)
sink()
##################