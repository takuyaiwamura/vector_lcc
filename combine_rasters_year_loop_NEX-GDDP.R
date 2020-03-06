## Combine the 8 tiles for all the years in a directory (loop)##
require(raster)
require(rasterVis)
require(rgdal)
require(geosphere)


TYPE <- "gen.pr"  #"gen.pr", "dev.m2.pr"
LCC.TYPE <-"pr"  #"lcc","ck","pr"

BOUNDARY <- "countries" # e.g. Africa, Europe, China, Mediterranan_wo_atlantic, North_America_wo_pacific, countries
BOUNDARY_DIR <- "../../../../GIS_dataset/countries/"       # Set GIS boundary directory
BOUNDARY.LOC <- paste(BOUNDARY_DIR,BOUNDARY,".shp",sep="")  # Africa., Europe., China., Mediterranan_wo_atlantic., North_America_wo_pacific.

SCENARIO <- "historical"  #"historical" "rcp45" "rcp85"
GCM <- "bcc-csm1-1"     #"bcc-csm1-1", "CCSM4","IPSL-CM5A-LR", "MIROC-ESM-CHEM"
RAIN <- "200"

RAS_DIR <- "out_rasters"
OUT_RAS_DIR <- "gcm_rasters"   # OUTPUT DIR (RASTER)
##OUT_FIG_DIR <- paste("gcm_figures",SCENARIO,sep="/")   # OUTPUT DIR (FIGURE)

## Read boundary map (for drawing maps)
boundary <-readOGR(dsn=BOUNDARY.LOC,layer=BOUNDARY)

## Create the list of rasters from the files to mach TYPE, REGION.NAME and YEAR
l <- list.files(RAS_DIR,full.names = TRUE)
l <- l[grep(TYPE,l)]
l <- l[grep(SCENARIO,l)]
l <- l[grep(GCM,l)]
l <- l[grep("Global-",l)]
l <- l[grep(paste("_",RAIN,".tif",sep=""),l)]
l <- l[grep("tif.",l,invert=TRUE)]  # remove anything with wierd file extensions e.g. tif.html
if (LCC.TYPE == "lcc"){
    l <- l[grep("_lcc_",l)]  # for LCC.TYPE = lcc
} else if (LCC.TYPE == "ck"){
    l <- l[grep("_ck_",l)]  # for LCC.TYPE = ck
} else {
    l <- l[grep("_ck",l,invert=TRUE)]  # for LCC.TYPE = pr
    l <- l[grep("_lcc",l,invert=TRUE)]  #for LCC.TYPE = pr
}

# Create year list from the file names
yr.l <- l[grep("Global-1",l)]
yr.l <- strsplit(yr.l,split="/")            # strip to file names
yr.l <- lapply(yr.l,'[[',length(yr.l[[1]]))
yr.l <- unlist(yr.l)
yr.l <- strsplit(yr.l,split="_")            # strip to years
yr.l <- lapply(yr.l,'[[',3)
yr.l <- unlist(yr.l)
yr.l <- as.numeric(yr.l)
yr.l <- yr.l[order(yr.l)]

print(paste(length(yr.l),l[1]))
                                        # Combine raster files (.tif) for each year

## MASK
mask.r <- raster("../global_average/continents2_id.tif")
mask.r <- mask.r / mask.r

for (i in 1:length(yr.l)) {
    print(yr.l[i])
    this.l <- l[grep(yr.l[i],l)]
    rlist <- lapply(this.l, raster)

    ## Project each raster to larger extent and the same origin.
    r.g <- extend(rlist[[1]],c(-180,180,-90,90))
    myProj <- function(x) {projectRaster(x, r.g, method='ngb')}
    rlist2 <- lapply(rlist,myProj)

    ## Merge all the rasters
    r.m <- do.call(raster::merge, rlist2)
    r.m.masked <- r.m * mask.r

    writeRaster(r.m.masked,paste(OUT_RAS_DIR,"/",TYPE,"_Global_",LCC.TYPE,"_",yr.l[i],"_",SCENARIO,"_",GCM,"_",RAIN,".tif",sep=""),'GTiff',overwrite=TRUE)

    ## ## Create a map with Option 1 (manual color schemes)
    ## breaks.gen <- c(0,0.01,5.01,10.01,15.01,20.01,25.01,40.01) # color breaks
    ## cols.gen <- c("lightgrey","navajowhite","khaki","yellow1","orange","tomato","red") # color breaks
    ## p1<-levelplot(r.m, col.regions=cols.gen, at=breaks.gen,margin=FALSE,main=paste("Number of life cycles in ",yr.l[i],sep=""),bg = "#A6CAE0")+layer(sp.lines(boundary, lwd=0.8, col='darkgray'))

    ## png(paste(OUT_FIG_DIR,"/",TYPE,"_Global_",yr.l[i],"_",SCENARIO,"_",GCM,".png",sep=""),width=1200,height=900)
    ## print(p1)
    ## dev.off()
}
