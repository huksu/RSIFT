require(jpeg)

gaussian2D <- function(x, y, sigma){
    return((1/(2*pi*sigma*sigma))*(exp(-(x*x + y*y)/(2*sigma*sigma))))
}

gaussian2DNonNorm <- function(x, y, sigma){
    # Normalization is not necessary now if it happens later
    return(exp(-(x*x + y*y)/(2*sigma*sigma)))
}

gaussian1DNonNorm <- function(x, sigma){
    # Normalization is not necessary now if it happens later
    return(exp(-(x*x)/(2*sigma*sigma)))
}

gaussianFilter2D <- function(sigma, radius = 0){
    radius <- ifelse(radius == 0, as.integer(sigma*3)+1, radius)
    filter <- matrix(data = 0, nrow = 2*radius-1, ncol = 2*radius-1)
    filter <- outer(1:nrow(filter),1:ncol(filter), FUN = Vectorize(function(r,c) gaussian2DNonNorm(radius-r,radius-c,sigma)))
    return(filter)
}

gaussianFilter1D <- function(sigma, radius = 0){
    radius <- ifelse(radius == 0, as.integer(sigma*3)+1, radius)
    filter <- c(-(radius-1):(radius-1))
    filter <- unlist(lapply(filter,FUN = function(x) gaussian1DNonNorm(x,radius)))
    return(filter)
}

simpleEdgeFilter2D <- function(){
    return(matrix(c(-1,-1,-1,-1,8,-1,-1,-1,-1), ncol = 3))
}

applyGaussFilter1D <- function(img, f2d){
    
    radius <- as.integer(nrow(f2d)/2)+1
    print(paste("Applying filter with radius:", radius))
    f <- f2d[,radius]
    
    # Col-wise
    newimg1 <- img
    newimg1 <- outer(1:nrow(newimg1),1:ncol(newimg1), FUN = Vectorize(function(i,j){
        imgcolsub <- j + c(1:length(f)) - radius
        imgcolsub <- ifelse(imgcolsub<1,1,imgcolsub)
        imgcolsub <- ifelse(imgcolsub>ncol(img),ncol(img),imgcolsub)
        return(sum(img[i,imgcolsub]*f))
    }))
    
    # Row-wise
    newimg2 <- newimg1
    newimg2 <- outer(1:nrow(newimg2),1:ncol(newimg2), FUN = Vectorize(function(i,j){
        imgrowsub <- i + c(1:length(f)) - radius
        imgrowsub <- ifelse(imgrowsub<1,1,imgrowsub)
        imgrowsub <- ifelse(imgrowsub>nrow(img),nrow(img),imgrowsub)
        return(sum(newimg1[imgrowsub,j]*f))
    }))
    
    return(imageNorm(newimg2))
}

applyGaussFilter2D <- function(img, f){
    # 2D filter leaves ugly borders on the shading.
    # Best to use 1D 2 pass filter.
    newimg <- img
    radius <- as.integer(nrow(f)/2)+1
    print(paste("Applying filter with radius:", radius))
    for(i in c(1:nrow(img))){
        for(j in c(1:ncol(img))){
            f_sub <- f[c(max(1,radius-i+1):min(nrow(f),nrow(img)-i+radius)),c(max(1,radius-j+1):min(ncol(f),ncol(img)-j+radius))]
            img_sub <- img[c(max(1,i-radius+1):min(nrow(img),i+radius-1)),c(max(1,j-radius+1):min(ncol(img),j+radius-1))]
            accumulation <- sum(as.vector(f_sub)*as.vector(img_sub))
            newimg[i,j] <- accumulation
        }
    }
    return(imageNorm(newimg))
}

imageNorm <- function(img){
    img <- img-min(img)
    img <- img/max(img)
    return(img)
}

imagePixClip <- function(img){
    img <- ifelse(img < 0, 0, img)
    img <- ifelse(img > 1, 1, img)
}

drawImg <- function(img){
    plot(c(0,1),c(0,1),t='n')
    rasterImage(img, 0,0,1,1)
}

grayImage <- function(colorImg){
    grayImg <- colorImg[,,1]+colorImg[,,2]+colorImg[,,3] # reduce to gray
    grayImg <- imageNorm(grayImg)
    return(grayImg)
}

downscaleImg <- function(img, level=2){
    halfheight <- floor(nrow(img)/level)
    halfwidth <- floor(ncol(img)/level)
    newimg <- matrix(0,nrow=halfheight,ncol=halfwidth)
    for(i in c(1:halfheight)){
        for(j in c(1:halfwidth)){
            newimg[i,j] <- img[(i-1)*level+1,(j-1)*level+1]
        }
    }
    return(newimg)
}

getK <- function(numScaleLevels){
    k <- 2^(1.0 / numScaleLevels)
}

getSigmas <- function(numScaleLevels, sd){
    sigmas <- rep(0,times=numScaleLevels)
    k <- getK(numScaleLevels)
    sigmas[1] <- sd
    sigmas[2] <- sd*sqrt((k^2)-1)
    for(i in c(3:(numScaleLevels+3))){
        sigmas[i] <- sigmas[i-1]*k
    }
    
    return(sigmas)
}

laplacianPyramid <- function(img, numOctaves = 4, numScaleLevels = 5, k = sqrt(2), sd = 1.6){
    
    sigmas <- getSigmas(numScaleLevels,sd)
    print("Sigma Values")
    print(sigmas)
    rootImg <- img
    
    lp <- list()
    
    for(i in c(1:numOctaves)){
        print(paste("Octave: ",i))
        octave <- array(rep(c(1:(numScaleLevels+3)), each = nrow(rootImg)*ncol(rootImg)), c(nrow(rootImg),ncol(rootImg),numScaleLevels+3))
        print(dim(octave))
        
        print(paste("Octave:", i, "Scale Level:", 1))
        octave[,,1] <- rootImg
        currentImg <- rootImg
        #drawImg(imageNorm(currentImg))
        
        for(j in c(2:(numScaleLevels+3))){
            print(paste("Octave:", i, "Scale Level:", j))
            f <- gaussianFilter2D(sigmas[j])
            nextImg <- applyGaussFilter1D(currentImg,f)
            #drawImg(imageNorm(nextImg))
            octave[,,j] <- nextImg
            currentImg <- nextImg
        }
        
        lp[[length(lp)+1]] <- octave
        
        # prep next loop
        rootImg <- downscaleImg(rootImg)
    }
    
    return(lp)
}

differenceOfGaussianPyramid <- function(lp, numOctaves, numScaleLevels){
    # The Difference-of-Gaussians pyramid is a list of 3d arrays
    # Each 3d array is a 2-D DoG image at different scale levels
    # Each 3d array in the list is the octave

    D <- list()
    for(i in c(1:numOctaves)){
        print(paste("Octave: ",i))
        rootImg <- lp[[i]][,,1]
        octave <- array(rep(c(1:(numScaleLevels+2)), each = nrow(rootImg)*ncol(rootImg)), c(nrow(rootImg),ncol(rootImg),numScaleLevels+2))
        print(dim(octave[,,]))
        
        for(j in c(1:(numScaleLevels+2))){
            print(paste("Octave:", i, "Scale Level:", j))
            diffImg <- imageNorm(lp[[i]][,,j+1] - lp[[i]][,,j])
            octave[,,j] <- imageNorm(diffImg)
        }
        
        D[[length(D)+1]] <- octave
    }
    
    return(D)
}

scaleSpacePeakDetection <- function(D, numOctaves, numScaleLevels, sd){
    peaks <- data.frame()
    neighbors <- do.call(expand.grid,list(x = -1:1, y = -1:1, z = -1:1))
    neighbors <- subset(neighbors, x != 0 | y != 0 | z != 0)
    for(o in c(1:numOctaves)){
        for(s in c(2:(numScaleLevels+1))){
            # for every scale (first and last notwithstanding), look for the extrema points
            print(paste("Octave:",o,"Scale:",s))
            for(r in c(2:(nrow(D[[o]][,,s])-1))){
                if(r%%10==0){
                    print(paste("Octave:",o,"Scale:",s,"Row:",r))
                }
                for(c in c(2:(ncol(D[[o]][,,s])-1))){
                    #print(paste("Col:",c))
                    # check every pixel's neighbors to see if this pixel is the maximum or minimum
                    greaterThanFail <- FALSE
                    lessThanFail <- FALSE
                    k <- 1
                    pixelValue <- D[[o]][r,c,s]
                    while((greaterThanFail == FALSE | lessThanFail == FALSE) & (k <= nrow(neighbors))){
                        neighborValue <- D[[o]][r+neighbors[k,]$x,c+neighbors[k,]$y,s+neighbors[k,]$z]
                        #print(paste("neighborValue:",neighborValue, "x:",r+neighbors[k,]$x, "x:",c+neighbors[k,]$y, "z:",s+neighbors[k,]$z))
                        greaterThanFail <- ifelse(pixelValue <= neighborValue,TRUE,greaterThanFail)
                        lessThanFail <- ifelse(pixelValue >= neighborValue,TRUE,greaterThanFail)
                        k <- k + 1
                    }
                    
                    if(greaterThanFail == FALSE | lessThanFail == FALSE){
                        # This point is an extrema
                        peaks <- rbind(peaks,data.frame(row = r, col = c, octave = o, interval = s))
                        print(paste("Extrema at r =", r, "c =",c))
                    }
                }
            }
        }
    }
    
    peaks$trow <- transformRCOPoint(peaks$row,peaks$octave)
    peaks$tcol <- transformRCOPoint(peaks$col,peaks$octave)
    peaks$scale <- peakPointScale(peaks$octave,peaks$interval,sd,numScaleLevels)
    peaks$octaveScale <- peakPointOctaveScale(peaks$interval,sd,numScaleLevels)
    
    return(peaks)
}

transformRCOPoint <- function(x,o){
    # we start with baselen at octave 1
    # when we downscale, we chew off pixels at the tail end
    # so we can simply multiply x without concern for middle points
    # 1 -> 1 -> 1
    # 2 -> 3 -> 5
    # 3 -> 5 -> 9
    # 4 -> 7 -> 13
    # 5 -> 9 -> 17
    # 6 -> 11 -> 21
    #... so on
    return((x*2^(o-1)) - (2^(o-1)) + 1)
}

rgbGrayImg <- function(grayImg){
    h <- dim(grayImg)[1]
    w <- dim(grayImg)[2]
    rgbGray <- array(0, dim = c(h,w,3))
    rgbGray[,,1] <- grayImg
    rgbGray[,,2] <- grayImg
    rgbGray[,,3] <- grayImg
    return(rgbGray)
}

setPixel <- function(img,x,y,r,g,b){
    if(x <= nrow(img) & y <= ncol(img) & x > 0 & y > 0){
        img[x,y,1] <- r
        img[x,y,2] <- g
        img[x,y,3] <- b
    }
    return(img)
}

setRGBGrayPeak <- function(rgbGray,r,c,o, drawCircles, radiusMult = 5){
    radius <- o*radiusMult
    
    # draw midpoint
    rgbGray <- setPixel(rgbGray,r,c,1,0,0)
    
    if(drawCircles){
        # draw circle
        for(x in c(-radius:radius)){
            y <- round(sqrt((radius * radius) - (x * x)))
            rgbGray <- setPixel(rgbGray,r+x,c+y,1,0,0)
            rgbGray <- setPixel(rgbGray,r+x,c-y,1,0,0)
        }
        for(y in c(-radius:radius)){
            x <- round(sqrt((radius * radius) - (y * y)))
            rgbGray <- setPixel(rgbGray,r+x,c+y,1,0,0)
            rgbGray <- setPixel(rgbGray,r-x,c+y,1,0,0)
        }    
    }
    
    return(rgbGray)
}

drawPeaks <- function(grayImg,peaks, drawCircles = TRUE){
    rgbGray <- rgbGrayImg(grayImg)
    for(i in c(1:nrow(peaks))){
        rgbGray <- setRGBGrayPeak(rgbGray,peaks[i,]$trow,peaks[i,]$tcol,peaks[i,]$octave,drawCircles)    
    }
    drawImg(rgbGray)
}


getRCSOffset <- function(D,r,c,s,o){
    # Get the derivative of D about r,c,s,o
    dx <- (D[[o]][r,c+1,s] - D[[o]][r,c-1,s]) / 2.0
    dy <- (D[[o]][r+1,c,s] - D[[o]][r-1,c,s]) / 2.0
    ds <- (D[[o]][r,c,s+1] - D[[o]][r,c,s-1]) / 2.0
    d <- matrix(c(dx,dy,ds),ncol=3)
    
    # Get the Hessian of D
    dxx <- D[[o]][r,c+1,s] + D[[o]][r,c-1,s] - 2 * D[[o]][r,c,s]
    dyy <- D[[o]][r+1,c,s] + D[[o]][r-1,c,s] - 2 * D[[o]][r,c,s]
    dss <- D[[o]][r,c,s+1] + D[[o]][r,c,s-1] - 2 * D[[o]][r,c,s]
    dxy <- (D[[o]][r+1,c+1,s] - D[[o]][r+1,c-1,s] - D[[o]][r-1,c+1,s] + D[[o]][r-1,c-1,s]) / 4.0
    dxs <- (D[[o]][r,c+1,s+1] - D[[o]][r,c-1,s+1] - D[[o]][r,c+1,s-1] + D[[o]][r,c-1,s-1]) / 4.0
    dys <- (D[[o]][r+1,c,s+1] - D[[o]][r-1,c,s+1] - D[[o]][r+1,c,s-1] + D[[o]][r-1,c,s-1]) / 4.0
    H <- matrix(c(dxx,dxy,dxs,dxy,dyy,dys,dxs,dys,dss), nrow=3, ncol=3)
    
    Hdet <- det(H)
    if(Hdet == 0){
        print("H determinant is 0.  Failing...")
        return(FALSE)
    }
    
    Hinv <- solve(H)
    
    rcsOffset <- round(as.vector(-Hinv %*% t(d)))
    names(rcsOffset)<-c("r","c","s")
    return(rcsOffset)
}

getContrast <- function(D,r,c,s,o,dr,dc,ds){
    # Get the derivative of D about r,c,s,o
    dx <- (D[[o]][r,c+1,s] - D[[o]][r,c-1,s]) / 2.0
    dy <- (D[[o]][r+1,c,s] - D[[o]][r-1,c,s]) / 2.0
    ds <- (D[[o]][r,c,s+1] - D[[o]][r,c,s-1]) / 2.0
    d <- matrix(c(dx,dy,ds),ncol=3)
    
    # Matrix form of xhat
    xhat <- matrix(c(dc,dr,dc), nrow=3)
    
    contrast <- D[[o]][r,c,s] + (d %*% xhat) * .5
    return(contrast)
}

peakLocalization <- function(peaks, D, numOctaves, numScaleLevels, sd, MAXSTEPS = 10, THRESHOLD = .03){
    # This is the part where most people get stuck on SIFT, so some explanation of what is happening:
    # We are going to wander around, starting from the current peak,
    # to find the point where the gradient of the pixel passes zero.
    # Imagine the pixel/image is represented by a function.  The gradient is the rate of 
    # change of that function.  Thinking of calculus, the point where the gradient
    # passes zero is the maximum or minimum of the original function.
    # Therefore, the derivative of the taylor series of D at 0 is used to approximate the gradient.
    # This approximation is the offset of the current pixel, because we are measuring the gradient from
    # the current pixel location.
    # We adjust x by the offset and repeat until the offset is not big enough to move the pixel anymore.
    # We also throw away localized peaks below a threshold.
    
    localizedPeaks <- data.frame()
    
    for(p in c(1:nrow(peaks))){
        print(paste("Peak: ",p))
        i <- 1
        success <- FALSE
        h <- dim(D[[peaks[p,]$octave]][,,peaks[p,]$interval])[1]
        w <- dim(D[[peaks[p,]$octave]][,,peaks[p,]$interval])[2]
        while(i <= MAXSTEPS){
            rcsOffset <- getRCSOffset(D,peaks[p,]$row,peaks[p,]$col,peaks[p,]$interval,peaks[p,]$octave)
            
            if(is.logical(rcsOffset)){
                print("RCSOffset Failed")
                break
            }
            
            if(sum(abs(rcsOffset)) == 0){
                # We didn't move, so finish
                print(paste("No movement at step:",i))
                success <- TRUE
                break
            }
            
            # Apply the offset
            peaks[p,]$row <- peaks[p,]$row + rcsOffset["r"]
            peaks[p,]$col <- peaks[p,]$col + rcsOffset["c"]
            peaks[p,]$interval <- peaks[p,]$interval + rcsOffset["s"]
            
            if(peaks[p,]$interval <= 1 | peaks[p,]$interval >= (numScaleLevels+2) |
                   peaks[p,]$row <= 1 | peaks[p,]$row >= h |
                   peaks[p,]$col <= 1 | peaks[p,]$col >= w ){
                # we went out of bounds
                print("Out of bounds...")
                break
            }
            
            i <- i + 1
        }
        
        if(i>MAXSTEPS){
            print("Too many steps...")
        }
        
        if(success == FALSE){
            # remove the point from the peak set
            print(paste("Fail on peak: ",p))
        } else {
            print(paste("Success on peak: ",p))
            contrast <- getContrast(D,peaks[p,]$row,peaks[p,]$col,peaks[p,]$interval,peaks[p,]$octave,rcsOffset["r"],rcsOffset["c"],rcsOffset["s"])
            if(1){ #TODO
                print(paste("Peak passed threshold test.  Contrast:",contrast))
                localizedPeaks <- rbind(localizedPeaks,peaks[p,])
            }
            else{
                print(paste("Peak failed threshold test.  Contrast:",contrast))
            }
            
        }
    }
    
    localizedPeaks$trow <- transformRCOPoint(localizedPeaks$row,localizedPeaks$octave)
    localizedPeaks$tcol <- transformRCOPoint(localizedPeaks$col,localizedPeaks$octave)
    localizedPeaks$scale <- peakPointScale(localizedPeaks$octave,localizedPeaks$interval,sd,numScaleLevels)
    localizedPeaks$octaveScale <- peakPointOctaveScale(localizedPeaks$interval,sd,numScaleLevels)
    return(localizedPeaks)
}

peakPointScale <- function(octave,interval,sd,numScaleLevels){
    return(sd * 2^(octave + interval/numScaleLevels))
}

peakPointOctaveScale <- function(interval,sd,numScaleLevels){
    return(sd * 2^(interval/numScaleLevels))
}

getStablePeaks <- function(peaks, D, numOctaves, numScaleLevels, radius = 10){
    stablePeaks <- data.frame()
    for(p in c(1:nrow(peaks))){
        print(paste("Peak: ",p))
        o <- peaks[p,]$octave
        r <- peaks[p,]$row
        c <- peaks[p,]$col
        s <- peaks[p,]$interval
        print(paste(o,r,c,s))
        
        dxx <- D[[o]][r,c+1,s] + D[[o]][r,c-1,s] - 2 * D[[o]][r,c,s]
        dyy <- D[[o]][r+1,c,s] + D[[o]][r-1,c,s] - 2 * D[[o]][r,c,s]
        dxy <- (D[[o]][r+1,c+1,s] - D[[o]][r+1,c-1,s] - D[[o]][r-1,c+1,s] + D[[o]][r-1,c-1,s]) / 4.0
        #H <- matrix(c(dxx,dxy,dxy,dyy), nrow=2, ncol=2)
        TraceH <- dxx + dyy
        DetH <- dxx*dyy - dxy^2
        test <- (TraceH^2)/DetH
        expected <- ((radius+1)^2)/radius 
        print(paste("Test: ",test,"Expected:",expected))
        if(DetH <= 0){
            print("Peak fail edge test. Determinant <= 0.")
        }
        else if(test <=  expected){
            print("Peak passed edge test.")
            stablePeaks<-rbind(stablePeaks,peaks[p,])
        }
        else {
            print("Peak failed edge test.  test > expected.")
        }
    }
    return(stablePeaks)    
}

assignOrientationAndMagnitude <- function(lp, peaks, sd, NUMBINS = 36, RADIUSFACTOR = 3, MAGTHRESHOLD = .8){ 
    
    radiusLevel <- sd * RADIUSFACTOR 
    orientedPeaks <- data.frame()
    peaks$orientation <- NA
    peaks$magnitude <- NA
    
    for(p in c(1:nrow(peaks))){
        print(paste("Orienting Peak: ",p))
        o <- peaks[p,]$octave
        r <- peaks[p,]$row
        c <- peaks[p,]$col
        s <- peaks[p,]$interval
        scale <- peaks[p,]$scale
        octaveScale <- peaks[p,]$octaveScale
        radius <- round(radiusLevel * octaveScale)
        print(paste("Pixel Radius:",radius))
        sigma <- RADIUSFACTOR * octaveScale
        
        # get the relevant picture from the gaussian pyramid
        img <- lp[[o]][,,s]
        h <- dim(img)[1]
        w <- dim(img)[2]
        
        hist <- rep(0,times = NUMBINS)
        
        for(i in c((-radius):radius)){
            for(j in c((-radius):radius)){
                # For each one of the surrounding points, look for the orientation and magnitude...
                # If you can calculate that point, then add it to the histogram...
                ri <- r + i
                ci <- c + i
                if(ri > 1 & ri < h & ci > 1 & ci < w){
                    # Calculate...
                    dx <- lp[[o]][ri,ci+1,s] - lp[[o]][ri,ci-1,s]
                    dy <- lp[[o]][ri-1,ci,s] - lp[[o]][ri+1,ci,s]
                    mag <- sqrt(dx*dx + dy*dy)
                    orient <- atan2(dy,dx)
                    
                    # Add to bin...
                    offsetDragFactor <- exp((-(i*i+j*j))/(2.0*sigma*sigma))
                    bin <- round(((orient + pi) / (2*pi)) * NUMBINS) + 1 # scale -PI,PI to 0,1
                    bin <- ifelse(bin > NUMBINS, 1, bin) # in case of 360 degrees
                    hist[bin] <- hist[bin] + mag * offsetDragFactor
                }
            }
        }
        
        # Find the max magnitude for each histogram and see if any of the other magnitudes
        # are above a threshold.  If they are, duplicate the peak.  Then set the orientations.
        maxVal <- 0
        for(bin in c(1:NUMBINS)){
            if(hist[bin] >= maxVal){
                maxVal <- hist[bin]
            }
        }
        
        print(paste("Max Magnitude:",maxVal))
        
        for(bin in c(1:NUMBINS)){
            if(hist[bin]/maxVal >= MAGTHRESHOLD){
                # Add the peak
                orientedPeaks <- rbind(orientedPeaks,peaks[p,])
                orientedPeaks[nrow(orientedPeaks),]$orientation <- ((2*pi*(bin-1))/NUMBINS) - pi
                orientedPeaks[nrow(orientedPeaks),]$magnitude <- hist[bin]
            }
        }
    }
    
    return(orientedPeaks)
}

generateSIFTDescriptors <- function(lp, orientedPeaks){
    SIFTDescriptors <- list()
    #for(p in c(1:nrow(orientedPeaks))){
    for(p in c(1:nrow(orientedPeaks))){
        
        print(paste("Describing Peak:",p))
        SUBREGIONW <- 4
        BINSPERSUBREGION <- 8
        RADIUSFACTOR <- 3
        
        o <- orientedPeaks[p,]$octave
        r <- orientedPeaks[p,]$row
        c <- orientedPeaks[p,]$col
        s <- orientedPeaks[p,]$interval
        scale <- orientedPeaks[p,]$scale
        octaveScale <- orientedPeaks[p,]$octaveScale
        sigma <- RADIUSFACTOR * octaveScale
        
        img <- lp[[o]][,,s]
        h <- dim(img)[1]
        w <- dim(img)[2]
        
        
        hist <-  array(0, dim = c(SUBREGIONW,SUBREGIONW,BINSPERSUBREGION))
        binSizeTheta <- BINSPERSUBREGION/(2*pi) # should be 25 degrees
        
        for(dr in c(-8:-1,1:8)){ # rows
            for(dc in c(-8:-1,1:8)){ # cols
                neighborr <- r + dr
                neighborc <- r + dc
                if(neighborr > 1 & neighborr < h & neighborc > 1 & neighborc < w){
                    # we can investigate this pixel
                    
                    if(neighborr < 0){
                        subregionr <- floor(dr/SUBREGIONW)+3
                    }
                    if(neighborr > 0){
                        subregionr <- ceiling(dr/SUBREGIONW)+2
                    }
                    if(neighborc < 0){
                        subregionc <- floor(dc/SUBREGIONW)+3
                    }
                    if(neighborc > 0){
                        subregionc <- ceiling(dc/SUBREGIONW)+2
                    }
                    
                    neighborr <- r + dr
                    neighborc <- r + dc
                    
                    # get the pixel orientation
                    dx <- lp[[o]][neighborr,neighborc+1,s] - lp[[o]][neighborr,neighborc-1,s]
                    dy <- lp[[o]][neighborr-1,neighborc,s] - lp[[o]][neighborr+1,neighborc,s]
                    mag <- sqrt(dx*dx + dy*dy)
                    orient <- atan2(dy,dx)
                    
                    # reorient the pixel to the feature orientation
                    orient <- orient + orientedPeaks[p,]$orientation
                    #print(paste("Orientation:",orient,orientedPeaks[p,]$orientation))
                    if(orient >= pi){
                        orient <- orient - 2*pi
                    }
                    
                    # Add to bin...
                    offsetDragFactor <- exp((-(dr*dr+dc*dc))/(2.0*sigma*sigma))
                    bin <- round(((orient + pi) / (2*pi)) * BINSPERSUBREGION) + 1 # scale -PI,PI to 0,1
                    bin <- ifelse(bin > BINSPERSUBREGION, 1, bin) # in case of 360 degrees
                    #print(paste("Bin:",bin, "subregionr:",subregionr,"subregionc:",subregionc))
                    #print(paste(mag,orient,bin,offsetDragFactor))
                    hist[subregionr,subregionc,bin] <- mag * offsetDragFactor
                    
                    
                }
            }
        }
        
        # Reorganize the bin weights into a single array
        descriptor <- list()
        for(x in c(1:SUBREGIONW)){
            for(y in c(1:SUBREGIONW)){
                for(z in c(1:BINSPERSUBREGION)){
                    descriptor[[length(descriptor)+1]] <- hist[x,y,z]
                }
            }
        }
        descriptor <- unlist(descriptor)
        maxmag <- max(descriptor)
        if(maxmag > 0){
            descriptor <- descriptor/maxmag # normalize to unit vector
        }
        SIFTDescriptors[[length(SIFTDescriptors)+1]] <- descriptor
    }
    
    return(SIFTDescriptors)
}

siftDataBuilder <- function(siftDescriptors, orientedPeaks){
    siftData <- list()
    for(p in c(1:nrow(orientedPeaks))){
        print(paste("Sift Point"))
        if(sum(siftDescriptors[[p]]) > 0){
            siftPoint <- list()
            siftPoint[[1]] <- orientedPeaks[p,]
            siftPoint[[2]] <- siftDescriptors[[p]]
            names(siftPoint) <- c("Peak","Descriptor")
            siftData[[length(siftData)+1]] <- siftPoint
        }
        else{
            print("Lack of magnitude in descriptor.")
        }
    }
    return(siftData)
}

SIFT <- function(imgFileName, numOctaves = 4, numScaleLevels = 5, sd = 1.6){
    ## Load the image
    myjpg <- readJPEG(imgFileName)
    grayImg <- grayImage(myjpg)
    
    ## Build the Difference-of-Gaussians pyramid
    lp <- laplacianPyramid(grayImg, numOctaves, numScaleLevels, sd)
    D <- differenceOfGaussianPyramid(lp, numOctaves, numScaleLevels)
    
    ## Scale space peak detection
    peaks <- scaleSpacePeakDetection(D, numOctaves, numScaleLevels, sd)
    
    ## Accurate keypoint localization
    localizedPeaks <- peakLocalization(peaks, D, numOctaves, numScaleLevels, sd, MAXSTEPS = 10)
    
    ## Eliminating Edge Responses
    stablePeaks <- getStablePeaks(localizedPeaks, D, numOctaves, numScaleLevels)
    
    ## Peak Orientation Assignment
    orientedPeaks <- assignOrientationAndMagnitude(lp, stablePeaks, sd)
    
    ## Generate the SIFT Descriptors (PHEW)
    drawPeaks(grayImg,orientedPeaks, drawCircles = FALSE)
    siftDescriptors <- generateSIFTDescriptors(lp, orientedPeaks)
    
    ## Combine the data into one set
    SIFTData <- siftDataBuilder(siftDescriptors,orientedPeaks)
    
    return(SIFTData)
}



SIFT.a <- SIFT("imgA.JPG")
SIFT.b <- SIFT("imgB.JPG")

length(SIFT.a)
length(SIFT.b)

