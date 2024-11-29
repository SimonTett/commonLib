# source("batch.filter_hourly.R")

# salloc -t 350 --mem=20000
# > Sys.setenv(RETICULATE_PYTHON = "/net/data/users/hsteptoe/miniconda3/envs/arrcc/bin/python")​
# > library(reticulate)​
# > use_python("/net/data/users/hsteptoe/miniconda3/envs/arrcc/bin/python")​
# > use_condaenv(condaenv="arrcc", required = T)​
# pre R
#     "module load scitools"
#     "conda activate sjb-iris-env"
# library(reticulate)
# iris <- import("iris")
# cube <- iris$load_cube("/scratch/hadsx/cpm/hourly/r001i1p00000/pr/TStest/r001i1p00000_20690201-20690230_pr.nc")
# v1 <- cube$data
# v1[300,1:654,1] <- NA
# cube$data <- v1
# iris$save(cube, '/scratch/hadsx/cpm/hourly/r001i1p00000/pr/TStest/r001i1p00000_20690201-20690230_pr.new.nc')


starttime <- Sys.time()

### test to see if in a BATCH environment and get pwd
cargs  <- commandArgs(trailingOnly = TRUE)
flist  <- cargs[1] # "/project/hires_rcm/UKCP18/cpm_rerun_output/TS3/hourly/r001i1p00000/pr/r001i1p00000_20690201-20690230_pr.nc"
stseas <- cargs[2] # 'djf'
# flist  <- "/scratch/hadsx/cpm/hourly/r001i1p00000/pr/TS1/r001i1p00000_19820601-19820630_pr.nc"
# stseas <- 'jja'

# flist  <- "/project/hires_rcm/UKCP18/cpm_rerun_output/TS3/hourly/r001i1p00000/pr/r001i1p00000_20680201-20680230_pr.nc"
# flist  <- "/scratch/hadsx/cpm/hourly/r001i1p00000/pr/TS3/r001i1p00000_20690201-20690230_pr.nc"
# flist  <- "/scratch/hadsx/cpm/hourly/r001i1p00000/pr/TStest/r001i1p00000_20690201-20690230_pr.nc"
# stout  <- sub('.nc','.new.nc',flist)
# stseas <- 'djf'

library(lattice)
library(fields)
library(Rfast)
# source(path("$RUTILS/simon_colours.R"))
# library(PCICt)
library(ncdf4)
source("/home/h03/hadsx/extremes/tawn/jordan/code/cpm/filter_fn.R")

cat("Loading", flist,cr)
# DIROUT0 <- '/scratch/hadsx/cpm/hourly/'
# dirbits <- strsplit(flist,'/')[[1]]
# DIROUT  <- paste(DIROUT0, dirbits[8], '/pr/', dirbits[6], '/', sep='')
# stout   <- paste(DIROUT, basename(flist),sep='')
# cat('Saving to',stout,cr)
cat('Saving to the input file',cr)

var_name  <- 'precipitation_flux' # 'unknown'
min_pr_th <- 10    # dont bother with fields whos max is below this value.  both files are in mm/h
rim       <- 54   # area to ignore around the rim
badb      <- 2 # no pixels to remove near bad data - big
bads      <- 1 # no pixels to remove near bad data - small
if(stseas=='djf') {
    badth   <- 180 # hbar or vbar over this is bad dta
    badthp  <- 12  # point over this is bad dta
    }
if(stseas=='jja') {
    badth   <- 160 #  hbar or vbar over this is bad dta
    badthp  <- 31  #  point over this is bad dta
}

nc1 <- nc_open(flist, write=TRUE)
v1  <- ncvar_get(nc1,var_name)
# if(nc1$var[[1]]$dim[[3]]$name!='time') stop("Time error")
# timeunits <- tail(unlist(strsplit(nc1$var[[1]]$dim[[3]]$units,'hours since ')),1)
# time      <- nc1$var[[1]]$dim[[3]]$vals
# calendar  <- nc1$var[[1]]$dim[[3]]$calendar
# pcictime  <- as.PCICt( time*(60*60), cal=calendar ,origin=timeunits)
# int_mon   <- as.integer(format(pcictime,"%m"))
# iseas     <- which(int_mon %in% imonyr[[stseas]] )
# v1        <- v1[,,iseas]
# > dim(v1) [1] 532 654 720

iroi.x    <- (1+rim):(dim(v1)[1]-rim)
iroi.y    <- (1+rim):(dim(v1)[2]-rim)

fmax <- 1:dim(v1)[3]
# for(i in 1:dim(v1)[3]) fmax[i]<-max(colMaxs(v1[iroi.x,iroi.y,i], value = TRUE, parallel = FALSE), na.rm=T)
# no Rfast version
for(i in 1:dim(v1)[3]) fmax[i]<-max(apply(v1[iroi.x,iroi.y,i],2,max,na.rm=TRUE), na.rm=T)

# readline("Stop")

# system.time(
for(i in 1:dim(v1)[3]){
# for(i in 1:20){

    # btime <- Sys.time()

    if(fmax[i]>min_pr_th) {

        # filter
        d1      <- filt_point_vbar_hbar_wnorm(v1[iroi.x,iroi.y,i], badth=badth, badthp=badthp, badb=badb, bads=bads)
        v1[iroi.x,iroi.y,i] <- d1
    }
    cat(i, ' ')
    # print(Sys.time() -btime)
}
# )

ncvar_put( nc1, var_name, v1, start=NA, count=NA, verbose=FALSE )

nc_close(nc1)

cat(cr,"Completed", basename(flist),cr)
print(Sys.time() -starttime)

#
