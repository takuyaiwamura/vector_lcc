#####  Phenology model  #####
##
## For A. aegypti
##
## Using NASA NEX-GDDP (0.25 degree/ 25 km); max temp, min temp, prcipitation
## NOTE: The NASA NEX-GDDP uses kg/m2/s as a unit for precipitation. 1 kg/m2/s = 86400 mm/day
##
## Mechanisms:
## Immature development GDD
## Egg hatching GDD
## Temperature sensitvite bloodmeal periods and gonotrophic cycle (GC)
## Heatkills and coldkills
## Precipitation constraints (annual)
## Cold constraints
##
## Features:
## Two running options: 1.Global tiles (REGION.ID 1-8). 2.Polygon (REGION.ID 0).
## Annual life cycles
## Monthly completion of develpment (egg hatching + immature development)
## Monthly completion of egg laying (bloodmeal + GC)
## Latitudinal average for monthly completions (only for global run)
##
## Outputs:
## Annual life cycles (generations)  --> Figures and Rasters
## Longest continous # of months above monthly development thresholds --> Figures and Rasters
## Monthly development completions (Egg hatching + Immature development) --> Figures
## Latitudinal average for monthly completions --> csv tables
##
##
##                                         by T. Iwamura  3/13/2019


## Load packages
require(rgdal)
require(ncdf4)
require(sp)
require(raster)
require(rasterVis)
require(lattice)
require(RColorBrewer)
require(geosphere)
require(maptools)
require(abind)
require(lubridate)
require(gridExtra)

###################
##  USER INPUTS  ##
###################

YEAR.FIRST <- 2050        # Year. Requires climate data for the YEAR and YEAR+1. (i.e. 2000 and 2001)
YEAR.LAST <- 2059        # Year. Requires climate data for the YEAR and YEAR+1. (i.e. 2000 and 2001)
SCE.ID <- 3         # 1:historical, 2:rcp45, 3:rcp85 (NOTE: Data range for 'historical': 1950-2005, for RCP 4.5/8.5: 2006-2100)
GCM <- "CCSM4" #"bcc-csm1-1", "CCSM4","IPSL-CM5A-LR", "MIROC-ESM-CHEM"

RAINFALL.TH <- 200

REGION.NAME <- "Global" # Used only for output file names
BOUNDARY <- "countries" # Name of shp file. (e.g. Africa, Europe, China, Mediterranan_wo_atlantic, North_America_wo_pacific, countries)

CLIMATE_DIR <- "../../../../GIS_dataset/NEX-GDDP/BCC/" # Set the data for NEX-GDDP data
BOUNDARY_DIR <- "../../../../GIS_dataset/countries/"       # Set GIS boundary directory

RAS_DIR <- "out_rasters/"   # OUTPUT DIR (RASTER)
TABLE_DIR <- "out_tables/"  # OUTPUT DIR (TABLE)
FIG_DIR <- "out_figures/"   # OUTPUT DIR (FIGURE)

## Below are autmatically created based on the parameters above
BOUNDARY.LOC <- paste(BOUNDARY_DIR,BOUNDARY,".shp",sep="")  # Africa., Europe., China., Mediterranan_wo_atlantic., North_America_wo_pacific.

if (SCE.ID == 1) {SCENARIO <- "historical"}
if (SCE.ID == 2) {SCENARIO <- "rcp45"}
if (SCE.ID == 3) {SCENARIO <- "rcp85"}

## NASA-NEX GDDP data selection is moved inside of the for loop (years).


###################
##  CUSTOM FUNC  ##
###################

## To convert an array to raster object
conv.a2r <- function (a,rncol,rnrow,ext,proj.str) {
    m <- matrix(a,ncol=rncol,nrow=rnrow)
    m <- t(m)
    m <- m[nrow(m):1,]
    r <- raster (m)
    r <- setExtent(r,ext)
    projection(r) <- proj.str
    return(r)
}

## To convert a raster object to an array
conv.r2a <- function (r,rncol,rnrow) {
    m <- raster::as.matrix(r)
    m <- m[nrow(m):1,]
    m <- t(m)
    a <- array(m)
    return(a)
}

## To convert an array to matrix (so that we can calc per lat/lon)
conv.a2m <- function (a,rncol,rnrow) {
    m <- matrix(a,ncol=rncol,nrow=rnrow)
    m <- t(m)
    m <- m[nrow(m):1,]
    return(m)
}

####################
##  MAIN PROGRAM  ##
####################
Sys.time()
for (y in YEAR.FIRST:YEAR.LAST){
    print(paste(y,":: ",Sys.time(),sep=""))
    YEAR <- y

    if (y == 2005){
        sce.2 <- "rcp45"
    } else {
        sce.2 <- SCENARIO
    }
    GDDP_MIN <- paste(CLIMATE_DIR,"tasmin/tasmin_day_BCSD_",SCENARIO,"_r1i1p1_",GCM,"_",YEAR,".nc",sep="")
    GDDP_MAX <- paste(CLIMATE_DIR,"tasmax/tasmax_day_BCSD_",SCENARIO,"_r1i1p1_",GCM,"_",YEAR,".nc",sep="")
    GDDP_MIN2 <- paste(CLIMATE_DIR,"tasmin/tasmin_day_BCSD_",sce.2,"_r1i1p1_",GCM,"_",YEAR+1,".nc",sep="")
    GDDP_MAX2 <- paste(CLIMATE_DIR,"tasmax/tasmax_day_BCSD_",sce.2,"_r1i1p1_",GCM,"_",YEAR+1,".nc",sep="")
    PRECIP <- paste(CLIMATE_DIR,"pr/pr_day_BCSD_",SCENARIO,"_r1i1p1_",GCM,"_",YEAR,".nc",sep="")


    ## OPEN netCDF4 data to collect info for creating arrays
    ## Get projection
    r <- raster(GDDP_MIN)
    proj.str <- proj4string(r)
    boundary <-readOGR(dsn=BOUNDARY.LOC,layer=BOUNDARY)

    ## GET data info for the nc data

    nc.min <- nc_open(GDDP_MIN)
    nc.max <- nc_open(GDDP_MAX)
    nc.min2 <- nc_open(GDDP_MIN2)
    nc.max2 <- nc_open(GDDP_MAX2)
    nc.pr <- nc_open(PRECIP)

    ## nc <- nc_open(GDDP_MIN)
    ## time <- ncvar_get(nc,"time")
    ## time.units <- ncatt_get(nc,"time","units")
    ## nt <- dim(time)
    ## name <- ncatt_get(nc,"tasmin","long_name")
    ## units <- ncatt_get(nc,"tasmin","units")
    ## fillValue <- ncatt_get(nc,"tasmin","_FillValue")
    ## missingValue <- ncatt_get(nc,"tasmin","missing_value")


    for (id in 1:8){
        if (id > 0) {
            out.name <- paste(REGION.NAME,"-", id, sep="")
            if (id == 1) { tile.ext = c(-179, -90, 0, 90) }
            if (id == 2) { tile.ext = c(-90, 0, 0, 90) }
            if (id == 3) { tile.ext = c(0, 90, 0, 90) }
            if (id == 4) { tile.ext = c(90, 179, 0, 90) }
            if (id == 5) { tile.ext = c(-179, -90, -90, 0) }
            if (id == 6) { tile.ext = c(-90, 0, -90, 0) }
            if (id == 7) { tile.ext = c(0, 90, -90, 0) }
            if (id == 8) { tile.ext = c(90, 179, -90, 0) }
        } else {tile.ext = c(0,0,0,0)}


        ## Set extent
        e <- tile.ext

        west <- e[1]
        east <- e[2]
        south <- e[3]
        north <- e[4]

        ## Set cell index for nc (regional)
        cell.south <- min(which(nc.min$dim$lat$vals < south + 0.2 & nc.min$dim$lat$vals > south - 0.2))
        cell.north <- min(which(nc.min$dim$lat$vals < north + 0.2 & nc.min$dim$lat$vals > north - 0.2))
        latIdx <- c(cell.south:cell.north) # update: use indexes instead of count

        west360 <- ifelse(west < 0, 360 + west,west)
        east360 <- ifelse(east < 0, 360 + east,east)

        cell.west <- min(which(nc.min$dim$lon$vals < west360 + 0.2 & nc.min$dim$lon$vals > west360 - 0.2)) + 1
        cell.east <- min(which(nc.min$dim$lon$vals < east360 + 0.2 & nc.min$dim$lon$vals > east360 - 0.2))
        if (cell.east > cell.west) {
            lonIdx <- c( cell.west:cell.east)
        } else {
            cell.360 <- which(nc.min$dim$lon$vals==max(nc.min$dim$lon$vals))
            lonIdx <- c(c(cell.west:cell.360),c(1:cell.east))
        }

        ##  Create array from the nc data (regional)
        ## MIN TEMP for the 1st year
        ##nc <- nc_open(GDDP_MIN)
        time <- ncvar_get(nc.min,"time")
        time.units <- ncatt_get(nc.min,"time","units")
        nt <- dim(time)
        name <- ncatt_get(nc.min,"tasmin","long_name")
        units <- ncatt_get(nc.min,"tasmin","units")
        fillValue <- ncatt_get(nc.min,"tasmin","_FillValue")
        missingValue <- ncatt_get(nc.min,"tasmin","missing_value")
        tasmin.array <- ncvar_get(nc.min,"tasmin")[lonIdx,latIdx,1:nt]
        tasmin.array[tasmin.array == fillValue$value] <- NA
        tasmin.array[tasmin.array == missingValue$value] <- NA
        ## MAX TEMP for the 1st year
        ##nc <- nc_open(GDDP_MAX)
        tasmax.array <- ncvar_get(nc.max,"tasmax")[lonIdx,latIdx,1:nt]
        tasmax.array[tasmax.array == fillValue$value] <- NA
        tasmax.array[tasmax.array == missingValue$value] <- NA
        ## MIN TEMP for the 2nd year
        ##nc <- nc_open(GDDP_MIN2)
        time2 <- ncvar_get(nc.min2,"time")
        time.units2 <- ncatt_get(nc.min2,"time","units")
        nt2 <- dim(time2)
        tasmin2.array <- ncvar_get(nc.min2,"tasmin")[lonIdx,latIdx,1:nt2]
        tasmin2.array[tasmin2.array == fillValue$value] <- NA
        tasmin2.array[tasmin2.array == missingValue$value] <- NA
        ## MAX TEMP for the 2nd year
        ##nc <- nc_open(GDDP_MAX2)
        tasmax2.array <- ncvar_get(nc.max2,"tasmax")[lonIdx,latIdx,1:nt2]
        tasmax2.array[tasmax2.array == fillValue$value] <- NA
        tasmax2.array[tasmax2.array == missingValue$value] <- NA

        ##  COMBINE MIN AND MAX AND FIX  ###
        ##Year 1
        tasmin.array <- tasmin.array - 273.15  ## convert K to C
        tasmax.array <- tasmax.array - 273.15  ## convert K to C
        tavg.array <- ( tasmin.array + tasmax.array ) / 2  ## Take average
        rm(tasmax.array, tasmin.array)
        ##Year 2
        tasmin2.array <- tasmin2.array - 273.15  ## convert K to C
        tasmax2.array <- tasmax2.array - 273.15  ## convert K to C
        tavg2.array <- ( tasmin2.array + tasmax2.array ) / 2  ## Take average
        rm(tasmax2.array, tasmin2.array)

        ## PRECIPITATION LAYER (regional)
        ##nc <- nc_open(PRECIP)
        pr.array <- ncvar_get(nc.pr,"pr")[lonIdx,latIdx,1:nt]
        time <- ncvar_get(nc.pr,"time")
        pr.time.units <- ncatt_get(nc.pr,"time","units")
        pr.nt <- dim(time)
        pr.name <- ncatt_get(nc.pr,"pr","long_name")
        pr.units <- ncatt_get(nc.pr,"pr","units")
        pr.fillValue <- ncatt_get(nc.pr,"pr","_FillValue")
        missingValue <- ncatt_get(nc.pr,"pr","missing_value")
        pr.array[pr.array == pr.fillValue$value] <- NA
        pr.array[pr.array == missingValue$value] <- NA
        pr.array <- pr.array * 86400 # convert from kg/m2/s to mm/d

        ## Convert to annual
        pr.y.array <- array(data=0,dim=c(nrow(pr.array[,,1]),ncol(pr.array[,,1])))
        for (i in 1:pr.nt) {
            pr.y.array <- pr.y.array + pr.array[,,i]
        }
        precip.y <- pr.y.array
        pr.y.array[pr.y.array<RAINFALL.TH] <- 0
        pr.y.array[pr.y.array>=RAINFALL.TH] <- 1
        pr.y.array <- array(pr.y.array)

        rm(pr.array)

        ## Convert nc ext (0 to 360) to boundary ext (-180 to 180)
        lon <- ncvar_get(nc.min,"lon")[lonIdx] ## subset with lat, lon
        lon <- ifelse(lon > 180, lon-360, lon) ## convert to minus longitude (-180 to +180)
        lat <- ncvar_get(nc.min,"lat")[latIdx] ## subset with lat, lon
        ext <- c(min(lon),max(lon),min(lat),max(lat))

        ## Create mask raster and mask array. Mask array is used to reduce comp time.
        dummy <- conv.a2r(array(tavg.array[,,1]),length(lat),length(lon),ext,proj.str)
        dummy <- setValues(dummy,NA)
        mask.r <- rasterize(boundary,dummy)
        mask.r[mask.r>0] <- 1
        mask.a <- conv.r2a(mask.r,length(lat),length(lon))  # masking array.

        ## SAVE WORKSPACE FOR BACKUP
        ## save.image("aegpti_precip.RData")

###############################
###                         ###
###      PHENOLOGY MODEL    ###
###                         ###
###############################

        ## Crete empty dataframe (model.df)
        ## NOTE: I changed the code for 'model.df' so that I don't need to have 40 sets of parameters. I added gen.n and changed the init for egg to 9999 (from NA).
        ## NOTE: I removed the egg.gen
        dl <- length(array(tavg.array[,,1]))
        model.df <- data.frame(array(data=NA,dim=dl),array(data=NA,dim=dl),  # $avg.t,$tmp.t
                               array(data=0,dim=dl),array(data=0,dim=dl),    # $gdd, $h.gdd
                               array(data=NA,dim=dl),                        # $hatch.days,$hatch
                               array(data=9999,dim=dl),                      # $imdev
                               array(data=4,dim=dl),array(data=8,dim=dl),    # $pre.blood,$gc
                               array(data=9999,dim=dl),                      # $egg
                               array(data=0,dim=dl),array(data=9999,dim=dl), # $cold.kill.days,$cold.kill
                               array(data=0,dim=dl),                         # $gen.n
                               array(data=0,dim=dl),array(data=0,dim=dl),    # $gdd.m, $h.gdd.m
                               array(data=0,dim=dl),                         # $egg.m
                               array(data=0,dim=dl)                          # $init.now
                               )

        colnames(model.df)<-c("avg.t","tmp.t",
                              "gdd","h.gdd",
                              "hatch",
                              "imdev",
                              "pre.blood", "gc",
                              "egg",
                              "cold.kill.days","cold.kill",
                              "gen.n",
                              "gdd.m","h.gdd.m",
                              "egg.m",
                              "init.now"
                              )

        record.out <- list()
        record.time <- 0

        ## GDD
        HATCH.BASE.TEMP <- 14.59   # Lower Temperature Threshold (from each referenced scenario)
        IMDEV.BASE.TEMP <- 11.78   # Lower Temperature Threshold (from each referenced scenario)

        HATCH.GDD <- 42.4          # Development Growing Degree Days (GDD) - completed immature development cycle (hatch to emerge)
        IMDEV.GDD <- 126.38        # Development Growing Degree Days (GDD) - completed immature development cycle (hatch to emerge)

        HEATKILL.THR <- 38.0       # HEAT.KILL conditions - 38. If the avg.t is above 40, then it kills development (convert DD to 0).
        COLDKILL.THR <- 0.0        # Cold kill threshold (if avg.t is below this, it kills development (convert DD to 0). If it's prolonged duration, they cannot colonize.
        COLDKILL.DAYS <- 152       # Cold days (define by COLDKILL.THR) required for cold kill condition completion

        ## Run the phenology model (loop: daily)
        this.month <- 1
        gdd.m.list <- list()
        h.gdd.m.list <- list()
        egg.m.list <- list()
        for (i in 1:nt){     # nt is the number of days in a year (nt is 365)
            ## Record data for monthly estimate
            model.df$egg.m <- ifelse(!model.df$egg == 9999, model.df$egg.m + (1 / model.df$egg), model.df$egg.m) # monthly development % for egg laying

            ## Update monthly generation number ($monthly.gen)...
            if (i == nt) {  ## If, i is the last day of the year,
                gdd.m.list[[this.month]] <- model.df$gdd.m      # gdd increase of this month
                h.gdd.m.list[[this.month]] <- model.df$h.gdd.m  # h.gdd increase of this month
                egg.m.list[[this.month]] <- model.df$egg.m      # egg % increase of this month

                print(paste("day:",i,"month:",this.month,"n egg.m:",length(model.df$egg.m)))

            } else if (month(as.Date(i-1,origin=as.Date(paste(YEAR,"01-01",sep="-")))) > this.month) { # If i finishes a month (i.e. i is the 1st day of a month),
                gdd.m.list[[this.month]] <- model.df$gdd.m      # gdd increase of this month
                h.gdd.m.list[[this.month]] <- model.df$h.gdd.m  # h.gdd  increase of this month
                egg.m.list[[this.month]] <- model.df$egg.m      # egg % increase of this month

                model.df$gdd.m <- ifelse(model.df$gdd.m>0, 0, 0)
                model.df$h.gdd.m <- ifelse(model.df$h.gdd.m>0, 0, 0)
                model.df$egg.m <- ifelse(model.df$egg.m>0, 0, 0)

                this.month <- this.month + 1

                print(paste("day:",i,"month:",this.month,"n egg.m:",length(model.df$egg.m)))
            }


            ## Prep
            model.df$avg.t <- array(tavg.array[,,i])       # this day's temperature
            model.df$tmp.t <- ifelse (model.df$avg.t >= 35, 35, model.df$avg.t) # if avg.t is above 35, there is no benefit

            ## Mask
            model.df$avg.t <- ifelse(is.na(mask.a), NA, model.df$avg.t) # if avg.t is above 35, there is no benefit
            model.df$tmp.t <- ifelse(is.na(mask.a), NA, model.df$tmp.t) # if avg.t is above 35, there is no benefit

            ## Hatching GDD
            model.df$h.gdd <- ifelse(model.df$avg.t >= HATCH.BASE.TEMP, model.df$h.gdd + model.df$tmp.t - HATCH.BASE.TEMP, model.df$h.gdd) # add GDD after hatching
            model.df$hatch <- ifelse((model.df$h.gdd >= HATCH.GDD) & is.na(model.df$hatch), i, model.df$hatch) # The day of hatching
            model.df$h.gdd <- ifelse(!is.na(model.df$hatch), 0, model.df$h.gdd)          # If hatched, reset $hatch.days to 0
            model.df$h.gdd <- ifelse (model.df$avg.t > HEATKILL.THR, 0, model.df$h.gdd) ## IF T is above HEATKILL.THR (e.g. 38C), reset GDD. (keep hatch condition 'on').

            ## Immature-development GDD
            ## ADD T  when T is above threshold.
            model.df$gdd <- ifelse(!is.na(model.df$hatch) & model.df$avg.t > IMDEV.BASE.TEMP, model.df$gdd + model.df$tmp.t - IMDEV.BASE.TEMP, model.df$gdd) # add GDD after hatching
            model.df$gdd <- ifelse (model.df$avg.t > HEATKILL.THR, 0, model.df$gdd) ## IF T is above HEATKILL.THR (e.g. 38C), reset GDD.
            model.df$hatch <- ifelse (model.df$avg.t > HEATKILL.THR, NA, model.df$hatch) ## IF T is above HEATKILL.THR (e.g. 38C), reset hatch
            model.df$gdd <- ifelse (model.df$avg.t < IMDEV.BASE.TEMP, 0, model.df$gdd) ## IF T is below IMDEV.BASE.TEMP (e.g. 12C), reset GDD.
            model.df$hatch <- ifelse (model.df$avg.t < IMDEV.BASE.TEMP, NA, model.df$hatch) ## IF T is below IMDEV.BASE.TEMP (e.g. 12C), reset hatch


            ## Monthly recroding
            model.df$h.gdd.m <- ifelse(is.na(model.df$hatch) & model.df$avg.t >= HATCH.BASE.TEMP, model.df$h.gdd.m + model.df$tmp.t - HATCH.BASE.TEMP, model.df$h.gdd.m)
            model.df$gdd.m <- ifelse(model.df$imdev == 9999 & !is.na(model.df$hatch) & model.df$avg.t > IMDEV.BASE.TEMP, model.df$gdd.m + model.df$tmp.t - IMDEV.BASE.TEMP, model.df$gdd.m)


            ## Record when the imature development is finished (=i), otherwise keep it as NA
            model.df$imdev <- ifelse((model.df$gdd > IMDEV.GDD) & model.df$imdev == 9999, i, model.df$imdev)

            ## Set pre blood conditions. Pre blood is set as 4.
            model.df$pre.blood <- ifelse (model.df$avg.t >= 16, 4, model.df$pre.blood) # if avg.t is above 16, set 4 days
            model.df$pre.blood <- ifelse (model.df$avg.t >= 20, 2, model.df$pre.blood) # if avg.t is above 20, set 2 days
            model.df$pre.blood <- ifelse (model.df$avg.t >= 26, 1, model.df$pre.blood) # if avg.t is above 26, set 1 day
            model.df$pre.blood <- ifelse (model.df$avg.t >= 35, 2, model.df$pre.blood) # if avg.t is above 35, set 2 days

            ## Set GC (Gonotrophic Cycle) conditions. GC is set as 8.
            model.df$gc <- ifelse (model.df$avg.t >= 20, 8, model.df$gc) # if avg.t is above 20, set 8 days
            model.df$gc <- ifelse (model.df$avg.t >= 26, 3, model.df$gc) # if avg.t is above 26, set 3 days
            model.df$gc <- ifelse (model.df$avg.t >= 30, 2, model.df$gc) # if avg.t is above 30, set 2 days
            model.df$gc <- ifelse (model.df$avg.t >= 35, 4, model.df$gc) # if avg.t is above 35, set 4 days

            ## Set the day when egg will be layed (pre.blood + gc) since $imdev. Otherwise keep it as 9999
            model.df$egg <- ifelse(model.df$imdev < 9999 & model.df$egg == 9999, model.df$pre.blood + model.df$gc, model.df$egg)

            ## Update:If day 'i' reaches the day of egg laying.. '$imdev+$egg',
            model.df$gen.n <- ifelse(i == model.df$imdev + model.df$egg, model.df$gen.n + 1, model.df$gen.n) # Increase generation n (life cycles)
            model.df$init.now <- ifelse(i == model.df$imdev + model.df$egg, 1, 0)                            # Activate initialization process

            ## Initialization
            model.df$h.gdd <- ifelse(model.df$init.now == 1, 0, model.df$h.gdd)           # Reset $h.gdd to 0
            model.df$gdd <- ifelse(model.df$init.now == 1, 0, model.df$gdd)               # Reset $gdd to 0
            model.df$hatch <- ifelse(model.df$init.now == 1, NA, model.df$hatch)          # Reset $hatch to NA
            model.df$imdev <- ifelse(model.df$init.now == 1, 9999, model.df$imdev)          # Reset $imdev to NA
            model.df$pre.blood <- ifelse(model.df$init.now == 1, 4, model.df$pre.blood)   # Reset $pre.blood to 4
            model.df$gc <- ifelse(model.df$init.now == 1, 8, model.df$gc)                 # Reset $gc to 8
            model.df$egg <- ifelse(model.df$init.now == 1, 9999, model.df$egg)            # Reset $egg to 9999
            model.df$init.now <- ifelse(model.df$init.now == 1, 0, 0)                     # Reset $init.now
        }

        ## Cold Kill (Egg death from cold kill during winter period of 2 years)
        tavg.2y <- abind(tavg.array,tavg2.array,along=3)  # need to read abind library
        nt2 <- dim(tavg.2y)[3]
        for (i in 180:540){   # nt is the number of days of the winter between YEAR and YEAR + 1
            model.df$avg.t <- array(tavg.2y[,,i])       # this day's temperature
            model.df$cold.kill.days <- ifelse((model.df$avg.t < COLDKILL.THR), model.df$cold.kill.days+1, 0)  # count days below avg temp.
            model.df$cold.kill <- ifelse((model.df$cold.kill.days == COLDKILL.DAYS) & (model.df$cold.kill == 9999), i, model.df$cold.kill) # length cold days
        }
        gen.ck <- ifelse(model.df$cold.kill==9999,model.df$gen.n,0)   # fix in 2019 - NA -> 0

        ## Remove the areas with too little precip (e.g. pr.y.array[pr.y.array>=200] <- 1)
        gen.pr <- ifelse(pr.y.array==1,gen.ck,0)        # fix in 2019 - NA -> 0

        ## Monthly development -> output to .csv table  ## REVISED
        dev.m <- list()
        dev.m.pr <- list()
        dev.m.avgs <- matrix(NA,nrow=length(lat),ncol=12)
        for (i in 1:12) {
            dev.m[[i]] <- (gdd.m.list[[i]]+h.gdd.m.list[[i]])/(HATCH.GDD+IMDEV.GDD)
            m.lcc <- conv.a2m(dev.m[[i]],length(lat),length(lon))
            dev.m.avgs[,i] <- apply(m.lcc,1,function (x) {mean(x,na.rm=TRUE)})
        }
        dev.m.avgs[dev.m.avgs=="NaN"] <- NA
        write.csv(dev.m.avgs, paste(TABLE_DIR,"dev_monthly_lcc_mean_",out.name,"_",YEAR,"_",SCENARIO,"_",GCM,"_",RAINFALL.TH,".csv",sep=""),row.names=FALSE)

        ## Monthly development -> output to .csv table  ## REVISED
        dev.m <- list()
        dev.m.ck <- list()
        dev.m.pr <- list()
        dev.m.avgs <- matrix(NA,nrow=length(lat),ncol=12)
        for (i in 1:12) {
            dev.m[[i]] <- (gdd.m.list[[i]]+h.gdd.m.list[[i]])/(HATCH.GDD+IMDEV.GDD)
            dev.m.ck[[i]] <- ifelse(model.df$cold.kill==9999,dev.m[[i]],0)     # fix in 2019 - NA -> 0
            dev.m.pr[[i]] <- ifelse(pr.y.array==1,dev.m.ck[[i]],0)     # fix in 2019 - NA -> 0

            m.pr <- conv.a2m(dev.m.pr[[i]],length(lat),length(lon))
            dev.m.avgs[,i] <- apply(m.pr,1,function (x) {mean(x,na.rm=TRUE)})
        }
        dev.m[dev.m=="NaN"] <- NA
        dev.m.ck[dev.m.ck=="NaN"] <- NA
        dev.m.pr[dev.m.pr=="NaN"] <- NA
        write.csv(dev.m, paste(TABLE_DIR,"dev_monthly_lcc_",out.name,"_",YEAR,"_",SCENARIO,"_",GCM,".csv",sep=""),row.names=FALSE)
        write.csv(dev.m.ck, paste(TABLE_DIR,"dev_monthly_ck_",out.name,"_",YEAR,"_",SCENARIO,"_",GCM,".csv",sep=""),row.names=FALSE)
        write.csv(dev.m.pr, paste(TABLE_DIR,"dev_monthly_pr_",out.name,"_",YEAR,"_",SCENARIO,"_",GCM,".csv",sep=""),row.names=FALSE)

        dev.m.avgs[dev.m.avgs=="NaN"] <- NA
        write.csv(dev.m.avgs, paste(TABLE_DIR,"dev_monthly_pr_mean_",out.name,"_",YEAR,"_",SCENARIO,"_",GCM,"_",RAINFALL.TH,".csv",sep=""),row.names=FALSE)

        ## ## Egg laying -> output to .csv table
        ## egg.m.pr <- list()
        ## egg.m.avgs <- matrix(NA,nrow=length(lat),ncol=12)
        ## for (i in 1:12) {
        ##     egg.m.ck[[i]] <- ifelse(model.df$cold.kill==9999,egg.m.list[[i]],NA)
        ##     egg.m.pr[[i]] <- ifelse(pr.y.array==1,egg.m.list[[i]],NA)
        ##     m.tmp <- conv.a2m(egg.m.pr[[i]],length(lat),length(lon))
        ##     egg.m.avgs[,i] <- apply(m.tmp,1,function (x) {mean(x,na.rm=TRUE)})
        ## }
        ## egg.m.avgs[dev.m.avgs=="NaN"] <- NA
        ## write.csv(egg.m.avgs, paste(TABLE_DIR,"egg_monthly_mean_",out.name,"_",YEAR,"_",SCENARIO,"_",GCM,".csv",sep=""),row.names=FALSE)


        ## Convert to raster ##
        ##  CONVERT arrays to master/raster (using 'raster' package)
        gen.n.r <- conv.a2r(model.df$gen.n,length(lat),length(lon),ext,proj.str)
        gen.ck.r <- conv.a2r(gen.ck,length(lat),length(lon),ext,proj.str)
        gen.pr.r <- conv.a2r(gen.pr,length(lat),length(lon),ext,proj.str)

        dev.m.pr.r = list()
        for (i in 1:12) {
            dev.m.pr.r[[i]] <- conv.a2r(dev.m.pr[[i]],length(lat),length(lon),ext,proj.str)
        }

        ## egg.m.pr.r = list()
        ## for (i in 1:12) {
        ##     egg.m.pr.r[[i]] <- conv.a2r(egg.m.pr[[i]],length(lat),length(lon),ext,proj.str)
        ## }


        ##  Raster outputs   ##
        writeRaster(gen.n.r, paste(RAS_DIR,"gen.pr_",out.name,"_lcc_",YEAR,"_",SCENARIO,"_",GCM,"_",RAINFALL.TH,".tif",sep=""),'GTiff',overwrite=TRUE)
        writeRaster(gen.ck.r, paste(RAS_DIR,"gen.pr_",out.name,"_ck_",YEAR,"_",SCENARIO,"_",GCM,"_",RAINFALL.TH,".tif",sep=""),'GTiff',overwrite=TRUE)
        writeRaster(gen.pr.r, paste(RAS_DIR,"gen.pr_",out.name,"_",YEAR,"_",SCENARIO,"_",GCM,"_",RAINFALL.TH,".tif",sep=""),'GTiff',overwrite=TRUE)

    }
    rm(nc.min,nc.min2,nc.max,nc.max2,nc.pr)
    gc()
}
Sys.time()

