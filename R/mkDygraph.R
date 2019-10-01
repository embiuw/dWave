#' \code{mkDygraph} Make interactive dygraph of data with (optional) detected lunges.
#' @param l.data List of data obtained from \code{run.d.wave}.
#' @return Creates an interactive plot where pre-determined lunges (if present), detected peaks (if present) and matched lunges & peaks.
#' @details Pre-determined lunges (if present) are indicated by red solid points.
#' Peaks detected using the 2D detector (if present) are represented by green solid points.
#' Matched lunges/peaks are enclosed in purple diamonds.
#' @family Lunge detector functions
#' @seealso \code{\link{run.d.wave}} to run \code{d.wave} on entire data record from an animal.
#' @author Martin Biuw
#' @examples
#' d <- read.depth(whalename='mn16_Jan25b')
#' t <- read.lunges(whalename='mn16_Jan25b')
#' df <- join.dstamp(d, t)
#' df$dive <- number.dive(df$depth)
#' df <- run.d.wave(df)
#' mkDygraph(df)
#' @export

mkDygraph <- function(l.data=tst) {
  require(dygraphs)
  require(xts)
  require(htmltools)

  if(class(l.data)!='list') {
    l.data <- list(df=l.data)
  }
  dyData <- data.frame(time=l.data$df$time, depth=l.data$df$depth)
  if(!is.na(match('lunge', names(l.data$df)))) {
    dyData$lunge <- dyData$depth
    dyData$lunge[which(!l.data$df$lunge)] <- NA
    include.lunges <- T
  } else {
    include.lunges <- F
  }
  if(!is.na(match('pks', names(l.data$df)))) {
    dyData$peak <- dyData$depth
    dyData$peak[which(!l.data$df$pks)] <- NA
    include.peaks <- T
  } else {
    include.peaks<- F
  }

  if(!is.na(match('matched', names(l.data$df)))) {
    dyData$matched <- dyData$depth
    dyData$matched[which(!l.data$df$matched)] <- NA
    include.matches <- T
  } else {
    include.matches <- F
  }

  if(!is.na(match('dive', names(l.data$df)))) {
    dAnn <- aggregate(l.data$df$time, list(l.data$df$dive),
                    function(x) mean(as.numeric(x)))
    names(dAnn) <- c('dive', 'time')
    dyData$dive <- rep(NA, nrow(dyData))
    ann.pos <- unlist(lapply(dAnn$time, function(x) {
      which.min(abs(x-as.numeric(dyData$time)))}))
    dyData$dive[ann.pos] <- min(dyData$depth)
    dAnn$time <- as.POSIXct(dAnn$time, origin='1970-01-01 00:00:00', tz='UTC')
    include.dive <- T
  } else {
    include.dive <- F
  }

  dts <- xts(dyData, order.by=dyData$time, tzone='UTC')
  yRange <- rev(range(dyData$depth))
  yRange <- yRange+(c(1, -1)*abs(0.04*diff(yRange)))


  dg <- dygraph(dts) %>%
    dyAxis('y', 'Depth (m)', valueRange=yRange) %>%
    dyOptions(drawGrid=F) %>%
    dyLegend(show='never')
  if(include.dive) {
    dySeries(dg, 'dive', drawPoints=T, pointSize=1)
  }
  if(include.matches) {
    dg <- dySeries(dg, 'matched', drawPoints=T, pointSize=8, pointShape='diamond', col='blue')
  }
  if(include.lunges) {
    dg <- dySeries(dg, 'lunge', drawPoints=T, pointSize=5, col='red')
  }
  if(include.peaks) {
    dg <- dySeries(dg, 'peak', drawPoints=T, pointSize=4, col='green')
  }
  dg <- dySeries(dg, 'depth', col='grey')
  dg <- dyRangeSelector(dg)
  if(include.dive) {
    for(i in 1:nrow(dAnn)) {
      eval(parse(text=paste('dg <- dg %>% dyAnnotation("', dAnn$time[i], '",
                            text=', dAnn$dive[i], ', series="dive", width=25)', sep='')))
    }
  }
    dg
}

