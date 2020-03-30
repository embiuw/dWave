#' \code{read.lunges} Read positions of lunges
#' @param filename Name of the Excel file containing lunge data
#' @param whalename Name of the sheet from which to read data
#' @return Returns a list with objects \code{start} (start time),
#' \code{freq} (sampling frequency) and \code{stamps} (index of stamp occurrences)
#' @details Note! R (or RStudio) must be run in 32-bit mode for this function to work.
#' @family Lunge detector functions
#' @seealso \code{\link{read.depth}} to read depth data,
#'   \code{\link{join.dstamp}} to merge lunge events and dive data,
#'   \code{\link{number.dive}} to add dive numeric ID column,
#' @author Martin Biuw
#' @examples
#' t <- read.lunges()
#' @importFrom XLConnect loadWorkbook readWorksheet
#' @export


read.lunges <- function(filename='./Data/Takashi data/Time stamp of Lunge for Norway.xlsx',
                        whalename='mn16_Jan19a') {

  require(XLConnect)
  xl <- loadWorkbook(filename)
  start <- readWorksheet(xl, whalename,
                         startRow=2, startCol=1, endRow=2, endCol=1,
                         header=F, colTypes='POSIXt')
  start <- as.POSIXct(as.numeric(start), origin='1970-01-01 00:00:00', tz='UTC')

  freq <- as.character(readWorksheet(xl, whalename,
                        startRow=2, startCol=2, endRow=2, endCol=2,
                        header=F, colTypes='character'))
  freq <- as.numeric(gsub(' Hz', '', freq))

  stamps <- as.character(readWorksheet(xl, whalename,
                                     startRow=2, startCol=3, endRow=0, endCol=3,
                                     header=F, colTypes='numeric'))
  stamps <- substr(stamps, 3, nchar(stamps)-1)
  stamps <- as.numeric(unlist(strsplit(stamps, ', ')))

  list(start=start, freq=freq, stamps=stamps)
}


#' \code{read.depth} Read depth data
#' @param whalename Name of the whale for which data are to be retrieved
#' @return Returns a list with objects \code{start} (start time),
#' \code{freq} (sampling frequency) and \code{d} (vector of depth data)
#' @details
#' @family Lunge detector functions
#' @seealso \code{\link{read.lunges}} to read lunge event timestamps,
#'   \code{\link{join.dstamp}} to merge lunge events and dive data,
#'   \code{\link{number.dive}} to add dive numeric ID column,
#' @author Martin Biuw
#' @examples
#' t <- read.depth()
#' @export

read.depth <- function(dirname='./Data/Takashi data/', whalename='mn16_Jan26', toff=1) {
  require(XLConnect)

  year <- unlist(strsplit(whalename, '_'))[1]
  year <- gsub('mn', '', year)
  year <- paste('20', year, sep='')

  if(!is.na(dirname)) {
    whaledir <- paste0(dirname, whalename)
    summaryfiles <- dir(dirname)[grep('Tagging', dir(dirname))]
  } else {
    whaledir <- paste(getwd(), whalename, sep='/')
    summaryfiles <- dir()[grep('Tagging', dir())]
  }

  filename <- paste(whaledir, '_Depth_ORG.TXT', sep='')
  which.summary <- summaryfiles[grep(year, summaryfiles)]
  if(is.na(dirname)) {
    xl <- loadWorkbook(which.summary)
    summ <- readWorksheet(xl, 1)
  } else {
    xl.file <- paste(dirname, which.summary, sep='/')
    xl.file <- gsub('//', '/', xl.file)
    xl <- loadWorkbook(xl.file)
    summ <- readWorksheet(xl, 1)
  }
  summ <- summ[match(tolower(whalename), tolower(summ$Data.ID)),]
  start <- as.POSIXct(as.numeric(summ$Record.start.time),
                                       origin='1970-01-01 00:00:00', tz='UTC')
  start <- start+(toff*3600)
  freq <- summ$sampling.inter.others

  d <- readLines(filename)
  if(d[1]=="Depth , Temp,  Speed,  Salinity and Admittance Data.") {
    d <- d[-c(1:grep('Depth', d)[2])]
    d <- unlist(lapply(d, function(x) {
      unlist(strsplit(x, ','))[1]
    }))
  }

  if(is.na(as.numeric(d[1]))) {
    while(is.na(as.numeric(d[1]))) {
      d <- d[-1]
    }
  }

  d <-as.numeric(d)
  if(abs(max(d))<abs(min(d))) {
    d <- -d
  }

  d <- d-min(d)
  d.offset <- sort(d)[which(cumsum(d)/max(cumsum(d))>0.5)[1]]
  d <- d-d.offset

  list(start=start, freq=freq, depth=d)
}


#' \code{join.dstamp} Merge depth and lunge data into data.frame
#' @param d Dive data, as created using \code{read.depth}
#' @param t Lunge data, as created using \code{read.lunge}
#' @return Returns a data frame with variables \code{time},
#' \code{depth} and \code{lunge}
#' @details
#' @family Lunge detector functions
#' @seealso \code{\link{read.lunges}} to read lunge event timestamps,
#'   \code{\link{read.depth}} to read dive data,
#'   \code{\link{number.dive}} to add dive numeric ID column,
#' @author Martin Biuw
#' @examples
#' t <- join.dstamp(d,t)
#' @export

join.dstamp <- function(d, t, trim=T) {
  d.int <- 1/d$freq
  t.int <- 1/t$freq
  d.df <- data.frame(time=d$start+(d.int*c(0:(length(d$depth)-1))), depth=d$depth)
  d.df$lunge <- rep(F, nrow(d.df))
  stamp.times <- t$start+((t.int*t$stamps)-t.int)
  stamp.match <- unlist(lapply(stamp.times, function(x) match(x, d.df$time)))
  d.df$lunge[stamp.match] <- T
  if(trim==T) {
    plot(-d.df$depth, type='l')
    loc <- round(locator(2)$x)
    if (loc[1]<1) loc[1] <- 1
    if(loc[2]>nrow(d.df)) loc[2] <- nrow(d.df)
    d.df <- d.df[c(loc[1]:loc[2]),]
  }
  plot(-depth~time, data=d.df, type='l')
  points(-depth~time, data=d.df[which(d.df$lunge),], col=2, xpd=NA)
  d.df
}


#' \code{dstamp.plot} Plot dive record with pre-determined lunges
#' @param df Data frame with data in dstamp format
#' @return Produces simple plot of dives and lunges
#' @details
#' @family Lunge detector functions
#' @seealso \code{\link{read.lunges}} to read lunge event timestamps,
#'   \code{\link{read.depth}} to read dive data,
#'   \code{\link{number.dive}} to add dive numeric ID column
#' @author Martin Biuw
#' @examples
#' dstamp.plot()
#' @export

dstamp.plot <- function(df) {
  plot(-depth~time, data=df, type='l')
  points(-depth~time, data=df[which(df$lunge),], col=2, xpd=NA)
}


#' \code{number.dive} Split data into distinct dives
#' @param data Depth variable in data frame created using \code{join.dstamp}
#' @param d.cut Initial cutoff (in meters) for dives to be detected.
#' @return Returns a data frame with variables \code{time},
#' \code{depth} and \code{lunge}
#' @details
#' @family Lunge detector functions
#' @seealso \code{\link{read.lunges}} to read lunge event timestamps,
#'   \code{\link{read.depth}} to read dive data,
#'   \code{\link{join.dstamp}} to add dive numeric ID column,
#' @author Martin Biuw
#' @examples
#' t <- number.dive(df)
#' @export

number.dive <- function(data=df$depth, d.cut=5) {
  data <- data-min(data)
  is.dive <- data>d.cut
  rl.dive <- rle(is.dive)
  rl.dive <- as.data.frame(do.call('cbind', rl.dive))
  rl.dive$values <- rl.dive$values==1
  rl.dive$starts <- c(1, cumsum(rl.dive$lengths)+1)[c(1:nrow(rl.dive))]
  which.dives <- which(rl.dive$values)
  dive.num <- rep(NA, length(data))
  for(i in 1:length(which.dives)) {
    this.dive <- rl.dive[which.dives[i],]
    inds <- c(this.dive$starts:(this.dive$starts+this.dive$lengths-1))
    diff.d <- data[inds[1]-1]-data[inds[1]]
    if(diff.d<=0) {
      while(diff.d<=0 & inds[1]>1) {
        inds <- c(inds[1]-1, inds)
        diff.d <- data[inds[1]-1]-data[inds[1]]
      }
      diff.d <- data[inds[1]]-data[inds[2]]
      if(diff.d==0) {
        while(diff.d==0) {
          inds <- inds[-1]
          diff.d <- data[inds[1]]-data[inds[2]]
        }
      }
    }
    diff.d <- data[inds[length(inds)]+1]-data[inds[length(inds)]]
    if(diff.d<=0) {
      while(diff.d<=0 & inds[length(inds)]<length(data)) {
        inds <- c(inds, inds[length(inds)]+1)
        diff.d <- data[inds[length(inds)]+1]-data[inds[length(inds)]]
      }
      diff.d <- data[inds[length(inds)]]-data[inds[length(inds)]-1]
      if(diff.d==0) {
        while(diff.d==0) {
          inds <- inds[-length(inds)]
          diff.d <- data[inds[length(inds)]]-data[inds[length(inds)]-1]
        }
      }
    }
    dive.num[inds] <- i
  }
  dive.num
}

#' \code{detect.lunge} Deprecated lunge detection function
#' @param
#' @return
#' @details
#' @family
#' @seealso
#' @author Martin Biuw
#' @examples
#' @export

detect.lunge <- function(data=dfc, method='resid',
                         window=4, square=F, q.pk=NA, a.pk=0.15, d.cut=0, pkdist=20,
                         plotting=F) {
  require(questionr)
  require(pracma)

  ddata <- diff(data$depth)
  p.var <- l.var <- power <- w.var <- w.mn <- vector('numeric')
  for(i in 1:length(ddata)) {
    wt <- dnorm(c(1:length(ddata)), i, window)
    wt <- wt/max(wt, na.rm=T)
    w.mn <- c(w.mn, wtd.mean(ddata, wt))
  }

  for(i in 1:length(ddata)) {
    if(!square) {
      wt <- dnorm(c(1:length(ddata)), i, window)
      wtl <- dnorm(c(1:length(ddata)), i, window*2)
      wt <- wt/max(wt, na.rm=T)
      wtl <- wtl/max(wtl, na.rm=T)
    } else {
      wt <- rep(0, length(ddata))
      which.wt <- c(max(c(1, i-window)):min(c(length(ddata),i+window)))
      wt[which.wt] <- 1
    }
    power <- c(power, sum(abs((ddata-w.mn)*wt)))
    p.var <- c(p.var, wtd.var((ddata-w.mn), wt, normwt=T))
    l.var <- c(l.var, wtd.var((ddata-w.mn), wtl, normwt=T))
  }

  if(method=='resid') {
    pk <- findpeaks(power, minpeakdistance=pkdist)
    if(is.na(q.pk)) {
      pk <- pk[which(pk[,1]>a.pk),]
      if(class(pk)!='matrix') pk <- matrix(pk, ncol=2)
    } else {
      pk <- pk[which(pk[,1]>quantile(pk[,1], q.pk)),]
      if(class(pk)!='matrix') pk <- matrix(pk, ncol=4)
    }
    if(any(pk[,2] %in% which(data$depth<d.cut))) {
      pk <- pk[-which(pk[,2] %in% which(data$depth<d.cut)),]
      if(class(pk)!='matrix') pk <- matrix(pk, ncol=4)
    }
  } else {
    pk <- findpeaks(p.var, minpeakdistance=pkdist)
    if(is.na(q.pk)) {
      pk <- pk[which(pk[,1]>a.pk),]
      if(class(pk)!='matrix') pk <- matrix(pk, ncol=4)
    } else {
      pk <- pk[which(pk[,1]>quantile(pk[,1], q.pk)),]
      if(class(pk)!='matrix') pk <- matrix(pk, ncol=4)
    }
    if(any(pk[,2] %in% which(data$depth<d.cut))) {
      pk <- pk[-which(pk[,2] %in% which(data$depth<d.cut)),]
    }
  }

  if(plotting) {
    par(mfrow=c(3,1), mar=c(1,4,1,1))
    plot(-data$depth, type='l', ylab='Depth (m)')
    points(which(data$lunge), -data$depth[which(data$lunge)], col=2, cex=2)
    plot(ddata, type='b', ylab=expression(paste(Delta, ' depth')))
    lines(w.mn, col=2)
    legend('bottomleft', lty=1, col=2,
           paste(window*2, '-sample gaussian weighted running mean'), bty='n')
    axis(2)
    box()
    if(method=='resid') {
      plot(power, type='n', ylab='Residual power')
      lines(abs(ddata-w.mn), type='l', col=3)
      lines(power, col=3)
    } else {
      plot(p.var, type='l', ylab='Smoothed local variance of residual')
      polygon(c(c(1:length(p.var)), rev(c(1:length(p.var)))),
              c(p.var, rep(0, length(p.var))), col=rgb(0,0,0,0.5))
    }
    points(pk[,2], pk[,1], col=2, cex=2)
    abline(v=pk[,2], col=2, lty=2, xpd=NA)
    par(mfrow=c(1,1))
  }
  if(method=='resid') {
    list(df=data.frame(w.mn=w.mn, resid=ddata-w.mn, power=power), peaks=pk)
  } else {
    list(df=data.frame(w.mn=w.mn, resid=ddata-w.mn, var=p.var), peaks=pk)
  }
}


#' \code{match.lunges} Lunge matching function
#' @param data Data frame with pre-determined lunges and detected peaks to be matched
#' @param window Period before and after detected peak within which pre-determined peaks are considered matched
#' @return list with 2 objects: \code{matched}, logical vector indicating which peaks were matched (or not) to a pre-determined lunge,
#' \code{confusion}, confusion matrix summarizing classification accuracy
#' @details NOTE! THis function is called by \code{run.d.wave}, so usually does not have to be run
#' @family Lunge detector functions
#' @seealso \code{\link{run.d.wave}} to run 2D detector on entire data record,
#' @author Martin Biuw
#' @examples
#' matched <- match.lunges(df)
#' @export

match.lunges <- function(data=df, window=15) {
  true.lunges <- which(data$lunge)
  if(length(true.lunges)>0) {
    matched <- unlist(lapply(true.lunges, function(x) {
      pk <- data$pks[c((x-window):(x+window))]
      any(pk)
    }))
    true.pos <- length(which(matched))
    false.neg <- length(which(!matched))
  } else {
    true.pos <- 0
    false.neg <- 0
  }

  obs.lunges <- which(data$pks)
  if(length(obs.lunges)>0) {
    matched <- unlist(lapply(obs.lunges, function(x) {
      pk <- data$lunge[c((x-window):(x+window))]
      any(pk)
    }))
    false.pos <- length(which(!matched))
  } else {
    false.pos <- length(obs.lunges)
  }

  no.lunges <- which(!data$lunge)
  unmatched <- unlist(lapply(no.lunges, function(x) {
    pk <- data$lunge[c((max(1, x-window)):(min(nrow(data), (x+window))))]
    any(pk)
  }))

  ## true.neg <- length(which(!unmatched))
  true.neg <- nrow(data)-sum(c(true.pos, false.pos, false.neg))

  conf.mat <- matrix(c(true.pos, false.pos, false.neg, true.neg), ncol=2)
  dimnames(conf.mat) <- list(Observed=c('True', 'False'), Predicted=c('True', 'False'))

  ## tpr <- conf.mat[1,1]/apply(conf.mat, 2, sum)[1]
  ## fnr <- conf.mat[2,1]/apply(conf.mat, 2, sum)[1]
  ## fpr <- conf.mat[1,2]/apply(conf.mat, 2, sum)[2]
  ## tnr <- conf.mat[2,2]/apply(conf.mat, 2, sum)[2]

  list(matched=matched, confusion=conf.mat)
}
