filt_point_vbar_hbar <- function(v1, badth=5, badthp=25, badb=4, bads=2 ) {

    # v1 2d image

    # depends(fields)
    # depends(raster)

    # badth    <- 5  # hbar or vbar over this is bad dta
    # badthp   <- 25 # point1 over this is bad dta
    # badb     <- 4  # no pixels to remove near bad data - big
    # bads     <- 2  # no pixels to remove near bad data - small

    vbar      <- matrix(1,5,3)
    vbar[,1]  <- -1
    vbar[,2]  <- 2
    vbar[,3]  <- -1
    vbar      <- vbar/length(vbar[])
    hbar      <- t(vbar)
    # for this definition of vbar and hbar somwhere around > +-5 is the right threshold
    pnt1      <- matrix(1/16,3,3)
    pnt1[2,2] <- 8/16
    mexpand   <- matrix(1,5,5)

    v2  <- t(v1)
    v2  <- raster::raster(v2[nrow(v2):1,]) # m <- t(m); raster(m[nrow(m):1,])
    v2v <- raster::focal(v2,vbar,pad=TRUE, na.rm=TRUE)
    v2h <- raster::focal(v2,hbar,pad=TRUE, na.rm=TRUE)
    v2p <- raster::focal(v2,pnt1,pad=TRUE, na.rm=TRUE)

    ### enlarge masked areas

    ### horizontal
    v2hm  <- raster::as.matrix(v2h)
    v2hm  <- v2hm[nrow(v2hm):1,]; v2hm <- t(v2hm)
    v2hm2 <- v2hm # NB overwriting v2hm
    for(d1 in seq_along(v2hm[,1])) {
        i5 <- which(v2hm2[d1,] < -badth | v2hm2[d1,] > badth)
        if(length(i5)>0) {
            for(d2 in i5) {
                x2 <- pmax((d1-badb),1):pmin((d1+badb),dim(v2hm)[1])
                y2 <- pmax((d2-bads),1):pmin((d2+bads),dim(v2hm)[2])
                v2hm[x2,y2] <- badth +1
            }
        }
    }
    # vertical
    v2vm  <- raster::as.matrix(v2v)
    v2vm  <- v2vm[nrow(v2vm):1,]; v2vm <- t(v2vm)
    v2vm2 <- v2vm # NB overwriting v2vm
    for(d1 in seq_along(v2vm[,1])) {
        i5 <- which(v2vm2[d1,] < -badth | v2vm2[d1,] > badth)
        if(length(i5)>0) {
            for(d2 in i5) {
                x2 <- pmax((d1-bads),1):pmin((d1+bads),dim(v2vm)[1])
                y2 <- pmax((d2-badb),1):pmin((d2+badb),dim(v2vm)[2])
                v2vm[x2,y2] <- badth +1
            }
        }
    }
    # point
    v2pm  <- raster::as.matrix(v2p)
    v2pm  <- v2pm[nrow(v2pm):1,]; v2pm <- t(v2pm)
    v2pm2 <- v2pm # NB overwriting v2pm
    for(d1 in seq_along(v2pm[,1])) {
        i5 <- which(v2pm2[d1,] < -badthp | v2pm2[d1,] > badthp)
        if(length(i5)>0) {
            for(d2 in i5) {
                x2 <- pmax((d1-1),1):pmin((d1+1),dim(v2pm)[1])
                y2 <- pmax((d2-1),1):pmin((d2+1),dim(v2pm)[2])
                v2pm[x2,y2] <- badthp +1
            }
        }
    }

    v1m      <- v1
    ixt      <- which(v2hm > badth | v2vm > badth | v2pm > badthp)
    v1m[ixt] <- NA

    v1m
}

### normalise wrt peak value in window
filt_point_vbar_hbar_wnorm <- function(v1, badth=240, badthp=46, badb=4, bads=2 ) {

    # v1 2d image

    # depends(fields)
    # depends(raster)

    # badth    <- 160 jja # 180 djf  # hbar or vbar over this is bad dta
    # badthp   <-  31 jja #  12 djf  # point1 over this is bad dta
    # badb     <- 2  # no pixels to remove near bad data - along main axis of bar
    # bads     <- 1  # no pixels to remove near bad data - perpendicular to main axix of bar

    ## define the convolution kernels
    vbar      <- matrix(1,5,3)
    vbar[,1]  <- -1
    vbar[,2]  <- 2
    vbar[,3]  <- -1
    vbar      <- vbar/length(vbar[])
    hbar      <- t(vbar)
    pnt1      <- matrix(-1/9,3,3)
    pnt1[2,2] <- 8/9
    # mexpand   <- matrix(1,5,5)

    ## convolve kernels onto data
    v2  <- t(v1)
    v2  <- raster::raster(v2[nrow(v2):1,]) # m <- t(m); raster(m[nrow(m):1,])
    v2v <- raster::focal(v2,vbar,pad=TRUE, na.rm=TRUE)
    v2h <- raster::focal(v2,hbar,pad=TRUE, na.rm=TRUE)
    v2p <- raster::focal(v2,pnt1,pad=TRUE, na.rm=TRUE)

    # ## normalise mk1 - found to be ineffective
    # v2hm  <- raster::as.matrix(v2h)
    # v2hm  <- v2hm[nrow(v2hm):1,]; v2hm <- t(v2hm)
    # iv10  <- which(v1>5)
    # x1    <- v2hm/v1
    # x1[-iv10] <- 0

    ## normalise mk2
    nhbar     <- hbar
    nhbar[]   <- 0
    nhbar[2,] <- 1

    ## find magnitude of rain in horisontal bars and scale by ^1/2
    v2hx   <- raster::focal(v2, w=nhbar, fun=max, pad=TRUE, na.rm=TRUE)
    v2hx[] <- v2hx[]^(1/2)
    v2hx[which(v2hx[]<1)] <- 1
    # normalise the convolution result
    nv2h <- 100*v2h/v2hx

    ## find magnitude of rain in verticle bars and scale by ^1/2
    v2vx <- raster::focal(v2, w=t(nhbar), fun=max, pad=TRUE, na.rm=TRUE)
    v2vx[] <- v2vx[]^(1/2)
    v2vx[which(v2vx[]<1)] <- 1
    # normalise the convolution result
    nv2v <- 100*v2v/v2vx

    ## not bothering to normalise the point convolution
    # v2px <- raster::focal(v2, w=matrix(1,3,3), fun=max, pad=TRUE, na.rm=TRUE)
    # v2px <- v2
    # v2px[] <- v2px[]^(1/6)
    # v2px[which(v2px[]<1)] <- 1
    # nv2p <- v2p/v2px
    nv2p <- v2p

    ### create masked fields and dialate by badb and bads

    ### horizontal
    v2hm  <- raster::as.matrix(nv2h)
    v2hm  <- v2hm[nrow(v2hm):1,]; v2hm <- t(v2hm)
    v2hm2 <- v2hm # NB overwriting v2hm
    for(d1 in seq_along(v2hm[,1])) {
        # i5 <- which(v2hm2[d1,] < -badth | v2hm2[d1,] > badth)
        i5 <- which(v2hm2[d1,] > badth)
        if(length(i5)>0) {
            for(d2 in i5) {
                x2 <- pmax((d1-badb),1):pmin((d1+badb),dim(v2hm)[1])
                y2 <- pmax((d2-bads),1):pmin((d2+bads),dim(v2hm)[2])
                v2hm[x2,y2] <- badth +1
            }
        }
    }
    # vertical
    v2vm  <- raster::as.matrix(nv2v)
    v2vm  <- v2vm[nrow(v2vm):1,]; v2vm <- t(v2vm)
    v2vm2 <- v2vm # NB overwriting v2vm
    for(d1 in seq_along(v2vm[,1])) {
        # i5 <- which(v2vm2[d1,] < -badth | v2vm2[d1,] > badth)
        i5 <- which(v2vm2[d1,] > badth)
        if(length(i5)>0) {
            for(d2 in i5) {
                x2 <- pmax((d1-bads),1):pmin((d1+bads),dim(v2vm)[1])
                y2 <- pmax((d2-badb),1):pmin((d2+badb),dim(v2vm)[2])
                v2vm[x2,y2] <- badth +1
            }
        }
    }
    # point
    v2pm  <- raster::as.matrix(nv2p)
    v2pm  <- v2pm[nrow(v2pm):1,]; v2pm <- t(v2pm)
    v2pm2 <- v2pm # NB overwriting v2pm
    for(d1 in seq_along(v2pm[,1])) {
        # i5 <- which(v2pm2[d1,] < -badthp | v2pm2[d1,] > badthp)
        i5 <- which(v2pm2[d1,] > badthp)
        if(length(i5)>0) {
            for(d2 in i5) {
                x2 <- pmax((d1-1),1):pmin((d1+1),dim(v2pm)[1])
                y2 <- pmax((d2-1),1):pmin((d2+1),dim(v2pm)[2])
                v2pm[x2,y2] <- badthp +1
            }
        }
    }

    # calculate combined mask for horizontal, vertical and point tests
    v1m      <- v1
    ixt      <- which(v2hm > badth | v2vm > badth | v2pm > badthp)
    v1m[ixt] <- NA

    return(v1m)
}


#######################################
find_vbar_hbar_wnorm <- function(v1, badth=240 ) {

    # v1 2d image

    # depends(fields)
    # depends(raster)

    # badth    <- 240 jja # 120 djf  # hbar or vbar over this is bad dta
    # badthp   <-  46 jja #   9 djf  # point1 over this is bad dta
    # badb     <- 4  # no pixels to remove near bad data - along main axis of bar
    # bads     <- 2  # no pixels to remove near bad data - perpendicular to main axix of bar

    vbar      <- matrix(1,5,3)
    vbar[,1]  <- -1
    vbar[,2]  <- 2
    vbar[,3]  <- -1
    vbar      <- vbar/length(vbar[])
    hbar      <- t(vbar)

    v2  <- t(v1)
    v2  <- raster::raster(v2[nrow(v2):1,]) # m <- t(m); raster(m[nrow(m):1,])
    v2v <- raster::focal(v2,vbar,pad=TRUE, na.rm=TRUE)
    v2h <- raster::focal(v2,hbar,pad=TRUE, na.rm=TRUE)

    # norm 1
    v2hm  <- raster::as.matrix(v2h)
    v2hm  <- v2hm[nrow(v2hm):1,]; v2hm <- t(v2hm)
    iv10  <- which(v1>5)
    x1    <- v2hm/v1
    x1[-iv10] <- 0

    # norm 2
    nhbar     <- hbar
    nhbar[]   <- 0
    nhbar[2,] <- 1

    v2hx   <- raster::focal(v2, w=nhbar, fun=max, pad=TRUE, na.rm=TRUE)
    # v2hx   <- v2
    v2hx[] <- v2hx[]^(1/2)
    v2hx[which(v2hx[]<1)] <- 1
    nv2h <- 100*v2h/v2hx

    v2vx <- raster::focal(v2, w=t(nhbar), fun=max, pad=TRUE, na.rm=TRUE)
    # v2vx <- v2
    v2vx[] <- v2vx[]^(1/2)
    v2vx[which(v2vx[]<1)] <- 1
    nv2v <- 100*v2v/v2vx

    # h.i      <- which(nv2h > badth)
    # v.i      <- which(nv2v > badth)

    # return(list(hval=nv2h, h.i=h.i, hval=nv2v, v.i=v.i))
    return(list(hval=abs(raster::as.matrix(nv2h)), vval=abs(raster::as.matrix(nv2v))))
}
