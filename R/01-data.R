library(mefa4)
library(QPAD)
library(maptools)
library(intrval)
library(raster)

## change this path according your recurring project location
od <- setwd("~/repos/recurring/offset")
load_BAM_QPAD(version = 3)
if (getBAMversion() != "3")
  stop("This script requires BAM version 3")
rlcc <- raster("./data/lcc.tif")
rtree <- raster("./data/tree.tif")
rtz <- raster("./data/utcoffset.tif")
rd1 <- raster("./data/seedgrow.tif")
crs <- proj4string(rtree)
source("functions.R")
setwd(od)


## load data for ON BBS

load("d:\\bam\\2021\\rof\\ON-BBS\\BBS.91.19.ON.stopsFeb2021.RData")
#load("d:\\bam\\2021\\rof\\ON-BBS\\BBS.ON.wide.RData")

m5$PKEY <- paste0(m5$SS, ":", m5$Year)
#compare_sets(m5$PKEY, spp_count$PKEY)

y1 <- Xtab(Abund ~ PKEY+Species_ID, m5)
x1 <- nonDuplicated(m5, PKEY, TRUE)[rownames(y1),]

head(x1$date)
head(x1$StartTime)
x1off <- make_x(
    dt=x1$date,
    tm=paste0(x1$StartTime.Hour, ":", x1$StartTime.Min),
    lon=x1$POINT_X,
    lat=x1$POINT_Y,
    dur=3,
    dis=Inf)
rownames(x1off) <- rownames(x1)
str(x)
summary(x1off)


s1 <- intersect(colnames(y1), getBAMspecieslist())

o1 <- matrix(0, nrow(x1off), length(s1))
rownames(o1) <- rownames(x1)
colnames(o1) <- s1

for (spp in s1) {
  cat(spp, "\n")
  flush.console()
  o1[,spp] <- make_off(spp, x1off)$offset
}
str(o1)
sum(is.na(o1))

## BAM V6
library(data.table)

