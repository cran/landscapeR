## ---- message=FALSE, warning=FALSE---------------------------------------
library(landscapeR)

library(raster)
m <- matrix(0, 33, 33)
r <- raster(m, xmn=0, xmx=10, ymn=0, ymx=10)

## ---- eval=FALSE---------------------------------------------------------
#  rr <- makePatch(r, size=500, rast=TRUE)
#  plot(rr)

## ---- warning=FALSE------------------------------------------------------
num <- 5
size <- 15
rr <- makeClass(r, num, size)
plot(rr)

## ---- warning=FALSE------------------------------------------------------
num <- 75
size <- 10
rr <- makeClass(r, num, size)
plot(rr)

## ------------------------------------------------------------------------
num <- 5
size <- c(1,5,10,20,50)
pts <- c(1, 33, 1089, 1057, 545)
rr <- makeClass(r, num, size, pts)
plot(rr)

## ---- eval=FALSE---------------------------------------------------------
#  rr <- makeClass(r, 1, size=500, rast=TRUE)
#  plot(rr)

## ------------------------------------------------------------------------
patchSize <- 500
newVal <- 3
centre <- 545
rr <- makeClass(r, 1, patchSize, centre, val=newVal, rast=TRUE)
plot(rr)

## ---- warning=FALSE------------------------------------------------------
rr <- makeClass(rr, 1, 100, bgr=newVal, rast=TRUE, val=5)
plot(rr)

## ------------------------------------------------------------------------
rr <- makeClass(r, 1, patchSize, centre, val=newVal, rast=TRUE)
plot(rr)
rr <- expandClass(rr, 3, 250)
plot(rr)

## ------------------------------------------------------------------------
m[,17] <- 1
r <- raster(m, xmn=0, xmx=10, ymn=0, ymx=10)
plot(r)
rr <- expandClass(r, 1, 200)
plot(rr)

## ------------------------------------------------------------------------
m[] <- 0
r <- raster(m, xmn=0, xmx=10, ymn=0, ymx=10)
rr = makeLine(r, size=50, val=2, convol=0.05, spt=545, rast=TRUE)
plot(rr)

