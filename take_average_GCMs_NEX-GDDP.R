## Ensemble (take average of) the output from different GCMs.

require(raster)
require(rasterVis)
require(rgdal)
require(geosphere)

TYPE <- "gen.pr"  #"gen.pr", "dev.m2.pr"

BOUNDARY <- "countries" # e.g. Africa, Europe, China, Mediterranan_wo_atlantic, North_America_wo_pacific, countries
BOUNDARY_DIR <- "../../../../GIS_dataset/countries/"       # Set GIS boundary directory
BOUNDARY.LOC <- paste(BOUNDARY_DIR,BOUNDARY,".shp",sep="")  # Africa., Europe., China., Mediterranan_wo_atlantic., North_America_wo_pacific.

SCENARIO <- "historical"  #"historical" "rcp45" "rcp85"

##RAS_DIR <- paste("gcm_rasters",SCENARIO,sep="/")   # INPUT DIR
RAS_DIR <- "gcm_rasters"

OUT_RAS_DIR <- "avg_rasters"   # OUTPUT DIR
##OUT_FIG_DIR <- paste("avg_figures",SCENARIO,TYPE,sep="/")   # OUTPUT DIR

RAIN <- 200
GCM.NUM <- 4

## Read boundary map (for drawing maps)
boundary <-readOGR(dsn=BOUNDARY.LOC)

## Create the list of rasters from the files to mach TYPE, REGION.NAME and YEAR
l <- list.files(RAS_DIR,full.names = TRUE)
l <- l[grep(SCENARIO,l)]
l <- l[grep(TYPE,l)]
l <- l[grep(paste(RAIN,".tif",sep=""),l)]
l <- l[grep("tif.",l,invert=TRUE)]  # remove anything with wierd file extensions e.g. tif.html


## Create year list from the file names
yr.l <- strsplit(l,split="/")
yr.l <- lapply(yr.l,'[[',length(yr.l[[1]]))  # take only filename
yr.l <- unlist(yr.l)
yr.l <- strsplit(yr.l,split="_")
yr.l <- lapply(yr.l,'[[',3)
yr.l <- unlist(yr.l)
yr.l <- as.numeric(yr.l)
yr.l <- yr.l[order(yr.l)]
yr.l <- unique(yr.l)

##print(paste(length(yr.l)))

# Combine raster files (.tif) for each year
for (i in 1:length(yr.l)) {
    print(yr.l[i])
    this.l <- l[grep(yr.l[i],l)]
    if (length(this.l) == GCM.NUM) {
        rlist <- lapply(this.l, raster)
        rstack <- stack(rlist)
        r.avg <- mean(rstack)

        writeRaster(r.avg,paste(OUT_RAS_DIR,"/",TYPE,"_Global_",yr.l[i],"_",SCENARIO,"_",RAIN,"_avg.tif",sep=""),'GTiff',overwrite=TRUE)

        ## ## Create a map with Option 1 (manual color schemes)
        ## if (TYPE == "gen.pr") {
        ##     breaks.gen <- c(0,0.01,5.01,10.01,15.01,20.01,25.01,40.01) # color breaks
        ##     fig.title <- paste("N of life cycles in ",yr.l[i],sep="")
        ## } else if (TYPE == "dev.m2.pr"){
        ##     breaks.gen <- c(0,0.001,1.01,3.01,5.01,8.01,10.01,12.01) # color breaks
        ##     fig.title <- paste("N of months with more than 2 cycles in ",yr.l[i],sep="")
        ## } else {
        ##     print("ERROR: INCORRECT TYPE")
        ## }
        ## cols.gen <- c("lightgrey","navajowhite","khaki","yellow1","orange","tomato","red") # color breaks
        ## p1<-levelplot(r.avg, col.regions=cols.gen, at=breaks.gen,margin=FALSE,main=fig.title,bg = "#A6CAE0")+layer(sp.lines(boundary, lwd=0.8, col='darkgray'))

        ## png(paste(OUT_FIG_DIR,"/",TYPE,"_Global_",yr.l[i],"_",SCENARIO,"_avg.png",sep=""),width=1200,height=900)
        ## print(p1)
        ## dev.off()
    } else {
        print(paste("Error:", yr.l[i], "layer number =", length(this.l)))
    }


}
