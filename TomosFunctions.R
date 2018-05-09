# tomosUtil functions
# includes some functions that are used often in my R scripts

# Tomo Eguchi
# 12 February 2014 Started


# Extracting posterior samples of deviance or any other variable from jags output:
extract.samples <- function(varname, zm){
  dev <- unlist(lapply(zm, FUN = function(x) x[, varname]))
  return(dev)
}

#` Computes DIC using Gelman's equation. Input is the summary of coda samples.
#` "deviance" should be one of the monitored parameters and "dic" module
#` needs to be loaded prior to defining model. load.module('dic') when
#` using rjags. 27 March 2018
#`
#` @param summary.zm summary of coda samples
#` @return a DIC value
#` @seealso
#` @export
#` @examples
#` jm <- jags.model(file = 'jags_model.txt',
#`                  data = bugs.data,
#`                  inits = inits.function,
#`                  n.chains = n.chains,
#`                  n.adapt = n.iter)
#` params <- c('c', 'theta', 'sigma', 'y', 'deviance')
#` zm <- coda.samples(jm,
#`                   variable.names = params,
#`                   n.iter = n.iter)
#` summary.zm <- summary(zm)
#` DIC(summary.zm)

DIC <- function(summary.zm){
	D <- summary.zm$statistics['deviance',]
	DIC <-  D['Mean'] + 0.5 * (D['SD'])^2
	names(DIC) <- NULL
	return(DIC)
}

#` Computes beta parameters from mean and variance
#`
beta.params <- function(m, v){
  a <- ((1 - m) * (m^2))/v - m
  a[a < 0] <- NA
  b <- a * (1 - m)/m
  ab <- list(alpha = a, beta = b)
  return(ab)
}

# gamma.params computes gamma parameters from mean and variance
gamma.params <- function(m, v){
  a <- (m^2)/v
  b <- m/v
  ab <- list(alpha = a, beta = b)
  return(ab)
}

# beta.stats computes mean and variance for beta distribution
beta.stats <- function(a, b){
  mean <- a/(a+b)
  var <- (a*b)/((a+b+1)*((a+b)^2))
  mode <- (1 - a)/(2 - a - b)
  stats <- list(mean = mean, var = var, mode = mode)
  return(stats)

}

# gamma.stats computes mean, variance and mode for gamma distribution
gamma.stats <- function(theta, kappa){
  mean <- theta / kappa
  var  <- theta / (kappa ^ 2)

  if (theta > 1){
    mode <- (theta-1)/kappa
  } else {
    mode <- NA
  }

  stats <- list(mean = mean, var = var, mode = mode)
  return(stats)
}


#getCoastLine
# A function to extract a coast line given longitude and latitude
# limits for the Pacific Ocean

# Data files - E and W Pacific are separated; got these from Rich Cosgrove

# Used to be in Matlab

# Tomo Eguchi
# 16 December 2016
getCoastLine <- function(filename, lon.limits, lat.limits){
  data1 <- read.table(filename, header = F, row.names = NULL, sep = ',')
  # fix longitudes to -180 to 180
  lon.limits[lon.limits > 180] <- lon.limits[lon.limits > 180] - 360
  out.df <- data1[data1[,1] >= min(lon.limits) &
                  data1[,1] <= max(lon.limits) &
                  data1[,2] >= min(lat.limits) &
                  data1[,2] <= max(lat.limits),]
  colnames(out.df) <- c('Longitude', 'Latitude')
  out.df$idx <- NA
  idx.rows <- 1:nrow(out.df)
  na.idx <- is.na(out.df$Longitude)  # T/F for NA
  dif.na.idx <- na.idx[2:length(na.idx)] - na.idx[1:(length(na.idx)-1)]
  idx.neg1 <- idx.rows[dif.na.idx == -1]  # beginning of segments
  idx.pos1 <- idx.rows[dif.na.idx == 1]   # end of segments

  for (k in 1:length(idx.neg1)) {
    out.df[idx.neg1[k]:idx.pos1[k], 'idx'] <- k
  }

  # change index to a factor variable
  out.df$idx <- as.factor(out.df$idx)
  out.df <- na.omit(out.df)
  # this splits all segments into separate dataframes
  out.list <- split(out.df, out.df$idx)

  for (k in 1:length(out.list)){
    line.1 <- out.list[[k]][1,]
    line.end <- out.list[[k]][nrow(out.list[[k]]),]
    if (line.1$Longitude != line.end$Longitude){
      tmp <- out.list[[k]]
      tmp <- rbind(c(max(tmp$Longitude), max(tmp$Latitude), k),
                 tmp, c(max(tmp$Longitude), max(tmp$Latitude), k))
      out.list[[k]] <- tmp
    }
  }
  return(out.list)

}

dirSelector <- function(){
  sysInfo <- Sys.info()

  # Change the directory structure according to the computer
  if (sysInfo[1] == 'Linux') {
    D00 <- '~/Documents/TomosFolders/'
    Rdir <- '~/Documents/R/'
    Dpub <- '~/Documents/TomosFolders/publications/'
  } else if (sysInfo[1] == 'Windows'){
    if (sysInfo[4] == 'SWC-TEGUCHI-D'){
      D00 <- "D:/TomosFolders/"
      Dpub <- "C:/Users/t_e/Documents/Publications/"
      Rdir <- "~/R/"
    }
    if (sysInfo[4] == 'SWC-TEGUCHI-L'){
      D00 <- "C:/Users/tomo.eguchi/Documents/TomosFolders/"
      Rdir <- "~/R/"
      Dpub <- "C:/Users/tomo.eguchi/Documents/Publications/"
    } else if (sysInfo[4] == 'SWC-TEGUCHI1-L'){
      D00 <- "D:/WorkStuff/"
      Rdir <- "~/R/"
      Dpub <- "D:/publications/"
    }

  }
  D <- list(Dtomo = D00, Rdir = Rdir, DPub = Dpub)
  return(D)
}

# find the equation of a line from two points
find.line <- function(xy1, xy2){
  #xy1 and xy2 are vectors of length 2, where the first
  # entry is an X coordinate and the second entry is
  # a Y coordinate.
  coef.b <- (xy1[2] - xy2[2])/(xy1[1] - xy2[1])
  coef.a <- xy1[2] - coef.b * xy1[1]
  coef <- list(slope = coef.b, intercept = coef.a)
  return(coef)
}

# finding an intersection between two lines given two coordinates per line
intersection2lines <- function(pts1, pts2){
  # pts1 and pts2 are 2x2 matrices, where each row of a matrix
  # is X and Y coordinates in that order.
  coef.line1 <- find.line(pts1[1,], pts1[2,])
  coef.line2 <- find.line(pts2[1,], pts2[2,])

  # create a matrix of coefficients. the second column is 1s
  # because y = a + bx is solved.
  A <- rbind(c(-coef.line1$slope, 1),
             c(-coef.line2$slope, 1))

  # check for the singularity of the coefficient matrix;
  # if it's singular, no intersection exists.
  QR <- qr(A)
  if (QR$rank == 1) {
    intersection <- NA
  } else {
    intersection <- solve(A) %*% c(coef.line1$intercept, coef.line2$intercept)
  }

  out <- c(X = intersection[1], Y = intersection[2])
  return(out)
}


# mm/dd/yyyy to DOY
mdY2DOY <- function(x) as.numeric(strftime(strptime(x, format="%m/%d/%Y"), "%j"))

# yyyy-mm-dd to DOY
YMD2DOY <- function(x) as.numeric(strftime(strptime(x, format="%Y-%m-%d"), "%j"))

# yyyy-mm-dd to Y
YMD2Y <- function(x) as.numeric(strftime(strptime(x, format="%Y-%m-%d"), "%Y"))

# yyyy-mm-dd to m
YMD2m <- function(x) as.numeric(strftime(strptime(x, format="%Y-%m-%d"), "%m"))

# yyyy-mm-dd to d
YMD2d <- function(x) as.numeric(strftime(strptime(x, format="%Y-%m-%d"), "%d"))

# mm/dd/yyyy to Year
mdY2Y <- function(x) as.numeric(strftime(strptime(x, format="%m/%d/%Y"), "%Y"))

# aa/bb/cc to c(aa, bb, cc)
splitByFwdSlash <- function(x) as.numeric(unlist(strsplit(as.character(x), '/')))

# yyyymmdd to Year
yyyymmdd2Y <- function(x) round(x/10000)

# yyyymmdd to Month
yyyymmdd2M <- function(x) round((x - 10000 * round(x/10000))/100)

# yyyymmdd to Day
yyyymmdd2D <- function(x) round(((x - 10000 * round(x/10000))/100 - round((x - 10000 * round(x/10000))/100)) * 100)

#yyyymmdd to DOY
yyyymmdd2DOY <- function(x){
  Y <- yyyymmdd2Y(x)
  m <- yyyymmdd2M(x)
  d <- yyyymmdd2D(x)
  YMD2DOY(paste0(Y, '-', m, '-', d))
}

mmddyy2date <- function(x){
  mm <- floor(x/10000)
  mm.char <- formatC(mm, width = 2, flag = '0')
  dd <- floor((x - mm*10000)/100)
  dd.char <- formatC(dd, width = 2, flag = '0')
  yy <- x - mm*10000 - dd*100
  Y <- ifelse(yy > 50, 1900 + yy, 2000 + yy)
  return(paste(Y, mm.char, dd.char, sep = "-"))
}

# hhmmss to hh:mm:ss  x is a character string or integer
hhmmss2hms <- function(x){
  library(lubridate)
  x.num <- as.numeric(x)
  hh <- floor(x.num/10000)
  mm <- floor((x.num - 10000 * hh)/100)
  ss <- round(((x.num - 10000 * hh)/100 - mm) * 100)
  return(paste0(formatC(hh, width = 2, flag = "0"), ':',
    formatC(mm, width = 2, flag = "0"), ':',
      formatC(ss, width = 2, flag = "0")))
}


# Multiple plot function
# from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# functions from http://rmflight.github.io/posts/2012/10/papersinRmd.html
incCount <- function(inObj, useName) {
    nObj <- length(inObj)
    useNum <- max(inObj) + 1
    inObj <- c(inObj, useNum)
    names(inObj)[nObj + 1] <- useName
    inObj
}

pasteLabel <- function(preText,
                       inObj,
                       objName,
                       postText = "",
                       insLink = TRUE) {
    objNum <- inObj[objName]
    useText <- paste(preText, objNum, ". ", postText, sep = "")
    if (insLink) {
        useText <- paste(useText, sep = "")
    }
    useText
}
## ###################

# finds figure/table number
labelNumber <- function(inObj, objName){
  objNum <- inObj[objName]
  return(objNum)
}

## convert nm to sm
nm2sm <- function(nm){
  sm <- 1.15078 * nm
  return(sm)
}

## convert sm to nm
sm2nm <- function(sm){
  nm <- sm / 1.15078
  return(nm)
}

## Convert nm to km
nm2km <- function(nm){
  km <- nm * 1.852
  return(km)
}

# convert km to nm
km2nm <- function(km){
  nm <- km / 1.852
  return(nm)
}

km2sm <- function(km){
  sm <- km/1.609
  return(sm)
}

sm2km <- function(sm){
  km <- sm * 1.609
  return(km)
}

deg2rad <- function(angleInDegrees){
  angleInRadians <- (pi/180) * angleInDegrees
  return(angleInRadians)
}

rad2deg <- function(angleInRadians){
  angleInDegrees <- angleInRadians * (180/pi)
  return(angleInDegrees)
}

ft2m <- function(distanceInFt){
  distanceInMeters <- distanceInFt * 0.3048
  return(distanceInMeters)
}

m2ft <- function(distanceInMeters){
  distanceInFt <- distanceInMeters/0.3048
  return(distanceInFt)
}

C2F <- function(x) x * 1.8 + 32

F2C <- function(x) (x - 32)/1.8

SE <- function(x) sd(x, na.rm = T)/sqrt(length(na.omit(x)))

# plotting gam output with ggplot - found it here:
# http://stackoverflow.com/questions/19735149/is-it-possible-to-plot-the-smooth-components-of-a-gam-fit-with-ggplot2
# removed ldply lines because they were not working AND I don't think
# they are necessary...
ggplot.model <- function(model,
                         type="conditional", res=FALSE,
                         col.line="black", col.point="black",
                         size.line=1, size.point=1) {
  require(visreg)
  require(plyr)
  plotdata <- visreg(model, type = type, plot = FALSE)
  smooths <- data.frame(Variable = plotdata$meta$x,
                        x = plotdata$fit[[plotdata$meta$x]],
                        smooth = plotdata$fit$visregFit,
                        lower = plotdata$fit$visregLwr,
                        upper = plotdata$fit$visregUpr)

  # smooths <- ldply(plotdata, function(part)
  #   data.frame(Variable = part$meta$x,
  #              x = part$fit[[part$meta$x]],
  #              smooth = part$fit$visregFit,
  #              lower = part$fit$visregLwr,
  #              upper = part$fit$visregUpr))
  #

  residuals <- data.frame(Variable = plotdata$meta$x,
                          x=plotdata$res[[plotdata$meta$x]],
                          y=plotdata$res$visregRes)

  # residuals <- ldply(plotdata, function(part)
  #   data.frame(Variable = part$meta$x,
  #              x=part$res[[part$meta$x]],
  #              y=part$res$visregRes))
  if (res)
    if (length(model$var.summary) > 1){
      p1 <- ggplot(smooths, aes(x, smooth)) +
            geom_line(col=col.line, size=size.line) +
            geom_line(aes(y=lower), linetype="dashed", col=col.line, size=size.line) +
            geom_line(aes(y=upper), linetype="dashed", col=col.line, size=size.line) +
            geom_point(data = residuals, aes(x, y), col=col.point, size=size.point) +
            facet_grid(. ~ Variable, scales = "free_x")
    } else {
      p1 <- ggplot(smooths, aes(x, smooth)) +
            geom_line(col=col.line, size=size.line) +
            geom_line(aes(y=lower), linetype="dashed", col=col.line, size=size.line) +
            geom_line(aes(y=upper), linetype="dashed", col=col.line, size=size.line) +
            geom_point(data = residuals, aes(x, y), col=col.point, size=size.point)
    }
  else
    if (length(model$var.summary) > 1){
      p1 <- ggplot(data = smooths, aes(x = x, y = smooth)) +
            geom_line(col=col.line, size=size.line) +
            geom_line(aes(y=lower), linetype="dashed", col=col.line, size=size.line) +
            geom_line(aes(y=upper), linetype="dashed", col=col.line, size=size.line) +
            facet_grid(. ~ Variable, scales = "free_x") +
            xlab(plotdata$meta$x)
    } else {
      p1 <- ggplot(data = smooths, aes(x = x, y = smooth)) +
            geom_line(col=col.line, size=size.line) +
            geom_line(aes(y=lower), linetype="dashed", col=col.line, size=size.line) +
            geom_line(aes(y=upper), linetype="dashed", col=col.line, size=size.line) +
            xlab(plotdata$meta$x)
    }

  return(p1)
}

inv.logit <- function(x) exp(x)/(1 + exp(x))
logit <- function(p) log(p/(1-p))

geo2math <- function(x){
  x <- x %% 360
  y <- 90 - x + 360 * ifelse(x>=270, 1, 0)
  return(y)
}

math2geo <- function(mdeg){
  gdeg <- 90 - mdeg
  gdeg[gdeg < 0] <- gdeg[gdeg < 0] + 360
  return(gdeg)

}

local2UTC <- function(local.time = Sys.time(), local.tz = Sys.timezone()){
  local <- as.POSIXct(local.time,
                      origin = '1970-01-01',
                      format = '%Y-%m-%d %H:%M:%S',
                      tz = local.tz)
  UTC <- lubridate::with_tz(local, tz = 'UTC')
  return(UTC)
}

UTC2local <- function(UTC.time, tz = Sys.timezone()){
  UTC.time <- as.POSIXct(UTC.time,
                           origin = '1970-01-01',
                           format = '%Y-%m-%d %H:%M:%S',
                           tz = 'UTC')
  local.time <- lubridate::with_tz(UTC.time, tz = tz)
  return(local.time)
}

remove.ext <- function(filename){
  parts <- unlist(strsplit(filename, '/'))
  filename <- parts[length(parts)]
  filename <- unlist(strsplit(filename, '\\.'))
  return(filename[1])
}

write.csv.rename <- function(x, file, quote = TRUE,
  eol = "\n", na = "NA", row.names = TRUE, fileEncoding = ""){

  if (file.exists(file)){
    file.creation.date <- format(file.info(file)$mtime, '%Y-%m-%d')
    if (length(grep(file, pattern = file.creation.date)) > 0){
      new.file.name <- paste0(unlist(strsplit(file, paste0(file.creation.date, '.csv'))),
                              file.creation.date, '.csv')
    } else {
      new.file.name <- paste0(unlist(strsplit(file, '.csv')),
                              '_', file.creation.date, '.csv')

    }

    file.rename(file, new.file.name)
  }

  write.csv(x, file = file, quote = quote,
      eol = eol, na = na, row.names = row.names,
      fileEncoding = fileEncoding)
}

replace.NA <- function(x, replaceWith = 0){
  # replaces NAs with a value - default is 0
  idx.NA <- is.na(x)
  x[idx.NA] <- replaceWith
  return(x)
}

# from https://stackoverflow.com/questions/19200841/consecutive-rolling-sums-in-a-vector-in-r
moving.cumsum <- function(x, n = 2){
  # computes cumulative sum over a span
  y <- rowSums(outer(1:(length(x)-n+1),
    1:n,
    FUN=function(i,j){x[(j - 1) + i]}))
  #y <- c(rep(NA, n-1), y)
  return(y)
}

plot.sst.SCB <- function(plot.date = "2018-01-01"){
  for (k in 1:length(plot.date)){
    if (as.Date(plot.date[k], format = "%Y-%m-%d") < as.Date("2017-10-01", format = "%Y-%m-%d")){
    URL <- paste0('http://coastwatch.pfeg.noaa.gov/erddap/griddap/jplG1SST.largePng?SST[(',
                plot.date[k], 'T00:00:00Z)][(31.005):(36.005)][(-121.995):(-116.995)]',
                '&.draw=surface&.vars=longitude%7Clatitude%7CSST',
                '&.colorBar=%7C%7C%7C12%7C26%7C&.bgColor=0xffccccff')

    download.file(URL,
                destfile = paste0('figures/jplG1SST_', plot.date[k], '.png'),
                mode='wb')

    } else {
      URL <- paste0("http://coastwatch.pfeg.noaa.gov/erddap/griddap/",
        "erdMWsstd14day.largePng?sst%5B(", plot.date[k], "T00:00:00Z)%5D%5B",
        "(0.0)%5D%5B(31.005):(36.005)%5D%5B(238.005):(243.005)%5D",
        "&.draw=surface&.vars=longitude%7Clatitude%7Csst&",
        ".colorBar=%7C%7C%7C%7C%7C&.bgColor=0xffccccff")
      download.file(URL,
                destfile = paste0('figures/MWsst14day_', plot.date[k], '.png'),
                mode='wb')

    }

  }
}

# converting three character month names to integer - may be in lubridate...
# added tolower() 27 March 2018
mmm2month <- function(x){
  switch(as.character(tolower(x)),
         "jan" = 1, "feb" = 2, "mar" = 3, "apr" = 4,
         "may" = 5, "jun" = 6, "jul" = 7, "aug" = 8,
         "sep" = 9, "oct" = 10, "nov" = 11, "dec" = 12, NA)
}
