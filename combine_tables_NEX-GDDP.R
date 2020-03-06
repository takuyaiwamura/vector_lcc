## Combine the 8 tiles for all the years in a directory (loop)##

##require(plotly)
##require(ggplot2)
##require(gridExtra)
##require(RColorBrewer)

TYPE <- "dev"

SCENARIO <- "historical"  #"historical" "rcp45" "rcp85"
TABLE.DIR <- "out_tables"
OUT.DIR <- "gcm_tables"
N.YEAR <- 5
GCMs <- c("bcc-csm1-1", "CCSM4","IPSL-CM5A-LR", "MIROC-ESM-CHEM")
RAIN <- 200
########### CREATE LAT AVERAGE ###########


for (i in 1:length(GCMs)){
    gcm <- GCMs[i]

    l <- list.files(TABLE.DIR,full.names = TRUE)
    l <- l[grep(TYPE,l)]
    l <- l[grep(gcm,l)]
    l <- l[grep("Global-",l)]
    l <- l[grep(SCENARIO,l)]
    l <- l[grep(paste(RAIN,".csv",sep=""),l)]

    ## Create year list from the file names
    yr.l <- l[grep("Global-1",l)]
    yr.l <- strsplit(yr.l,split="/")            # strip to file names
    yr.l <- lapply(yr.l,'[[',length(yr.l[[1]]))
    yr.l <- unlist(yr.l)
    yr.l <- strsplit(yr.l,split="_")            # strip to years
    yr.l <- lapply(yr.l,'[[',6)
    yr.l <- unlist(yr.l)
    yr.l <- as.numeric(yr.l)
    yr.l <- yr.l[order(yr.l)]

    for (y in 1:length(yr.l)){
        print(paste(SCENARIO,gcm,yr.l[y]))
        ## Create the list of rasters from the files to mach TYPE, REGION.NAME and YEAR
        this.files <- l
        this.files <- this.files[grep(yr.l[y],this.files)]

        n.hemi <- NULL
        for (j in 1:4) {
            this.dat <- as.matrix(read.csv(this.files[j]))
            this.dat <- as.vector(this.dat)
            if(is.null(n.hemi)){
                n.hemi <- this.dat
            }else{
                n.hemi <- rbind(n.hemi,this.dat)
            }
        }
        n.hemi.mean <- colMeans(n.hemi,na.rm=TRUE)
        n.hemi.mean <- matrix(n.hemi.mean,ncol=12)

        s.hemi <- NULL
        for (j in 6:8) {
            this.dat <- as.matrix(read.csv(this.files[j]))
            this.dat <- as.vector(this.dat)
            if(is.null(s.hemi)){
                s.hemi <- this.dat
            }else{
                s.hemi <- rbind(s.hemi,this.dat)
            }
        }
        s.hemi.mean <- colMeans(s.hemi,na.rm=TRUE)
        s.hemi.mean <- matrix(s.hemi.mean,ncol=12)

        global.mean <- rbind(n.hemi.mean,s.hemi.mean)
        colnames(global.mean) <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
        write.csv(global.mean,paste(OUT.DIR,"/dev_monthly_pr_mean_Global_",yr.l[y],"_",SCENARIO,"_",gcm,"_",RAIN,".csv",sep=""),row.names=FALSE)
    }
}


