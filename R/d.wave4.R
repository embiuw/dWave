#' \code{d.wave} 2D lunge detector for individual dive, with validation
#' @param data Data frame containing dive info (e.g. created using \code{read.depth}, \code{read.lunges},
#' \code{join.dstamp} and \code{number.dive}).
#' @param dive Sequence number of dive to be analysed.
#' @param presmooth.window Initial window of smooth for gaussian running mean filter
#' @param pk.per Peak wavelet period for weighted power analysis
#' @param per.rng Range of all periods to be considered in weighted power analysis
#' @param min.pkdist Minimum period between consecutive peaks
#' @param a.pk Minimum amplitude for peak to be considered significant
#' @param pk.env Proportion of a.pk, used to set the duration of what is considered a peak
#' @param return.wavelet Logical, whether wavelet object should be returned or only the detected peaks
#' @param outliers Logical, either flag outliers or delete them (flag by default, change to 'omit' to omit)
#' @param edges Logical, either flag edge peaks or delete them (flag by default, change to 'omit' to omit)
#' @param plotting Whether to plot the results or not
#' @return Returns a list with objects \code{parameters}, containing all parameters used when running detector,
#' \code{wavelet} containing the wavelet output (if \code{return.wavelet} was set to TRUE) and \code{peaks} data frame containing the detected peaks
#' @details
#' @family Lunge detector functions
#' @seealso \code{\link{read.lunges}} to read lunge event timestamps,
#'   \code{\link{read.depth}} to read depth time series,
#'   \code{\link{join.dstamp}} to merge lunge events and dive data,
#'   \code{\link{number.dive}} to add dive numeric ID column,
#'   \code{\link{run.d.wave}} to run \code{d.wave} on entire data record from an animal.
#' @author Martin Biuw
#' @examples
#' d <- read.depth(whalename='mn16_Jan25b')
#' t <- read.lunges(whalename='mn16_Jan25b')
#' df <- join.dstamp(d, t)
#' df$dive <- number.dive(df$depth)
#' d103 <- d.wave(dive=103)
#' @export


d.wave <- function(data=df, dive=96, presmooth.window=4,
                    pk.per=16, per.rng=20,
                    min.pkdist=40, a.pk=0.05, pk.env=0.85, return.wavelet=F,
                    outliers='flag', edges='flag', plotting=T) {

  require(WaveletComp)
  require(questionr)
  require(pracma)
  require(outliers)

  options(warn=-1)

  parameters <- list(dive=dive, presmooth.window=presmooth.window,
                     pk.per=as.numeric(pk.per), per.rng=as.numeric(per.rng),
                     min.pkdist=min.pkdist, a.pk=a.pk)

  dat <- data[which(data$dive==dive),]
  dat$d.depth <- c(diff(dat$depth), NA)
  dat$w.mn <- rep(NA, nrow(dat))

  ## Gaussian smooth on d.depth:
  for(i in 1:nrow(dat)) {
    wt <- dnorm(c(1:nrow(dat)), i, presmooth.window)
    wt <- wt/max(wt)
    dat$w.mn[i] <- wtd.mean(dat$d.depth, wt, na.rm=T)
  }

  ## Residual after smooth:
  dat$r.depth <- dat$d.depth-dat$w.mn

  ## Run continuous wavelet analysis on residual:
  dw <- analyze.wavelet(dat[which(!is.na(dat$d.depth)),],
                        make.pval=F, 'r.depth', loess.span=0,
                        verbose=F)

  pers <- which.min(abs(dw$Period-pk.per))
  pers <- dw$Period[pers+(c(0,-per.rng, per.rng))]

  per.rrng <- range(pers)
  pkp <- dw$Period[which.min(abs(pk.per-dw$Period))]

  dr <- reconstruct(dw, plot.rec =F, rescale=F, siglvl=0.1,
                    sel.lower=per.rrng[1], sel.upper=per.rrng[1],
                    only.coi = T, verbose=F)

  r.sig <- dr$series$r.depth.r

  periods <- dw$Period

  lper <- log(periods)
  l.pk <- log(pkp)
  l.mn <- log(per.rrng[1])
  l.mx <- log(per.rrng[2])

  ## Set weighting sd to ensure weights drop to
  ## small value at mean and max periods
  sd <- optimize(function(sd) {
    min(abs(pnorm(l.mn, l.pk, sd)-0.001))
  }, interval=c(0.001, 1))$minimum

  ## Weighted power average, centred on peak period of 16s
  w.pow <- apply(dw$Power, 2, function(x)
    wtd.mean(x, dnorm(log(dw$Period), l.pk, sd), normwt=T))

  fp <- findpeaks(w.pow, minpeakdistance=min.pkdist)
  fp <- fp[order(fp[,2]),]
  if(class(fp)!='matrix') {
    fp <- matrix(fp, ncol=4)
  }
  fp <- data.frame(fp)

  names(fp) <- c('power', 'peak', 'start', 'end')

  fp <- fp[which(fp$power>=a.pk),]

  for(i in 1:nrow(fp)) {
    fp$start[i] <- fp$start[i] +
      which(w.pow[c(fp$start[i]:fp$peak[i])]>pk.env*fp$power[i])[1]-1
    fp$end[i] <- fp$peak[i] +
      tail(which(w.pow[c(fp$peak[i]:fp$end[i])]>pk.env*fp$power[i]), 1)-1
  }

  fp$outliers <- fp$edges <- rep(F, nrow(fp))

  ## Detect (& possibly omit) edge peaks:
  edge <- which(fp[,2]<=pk.per | fp[,2]>=(length(w.pow)-pk.per))
  if(length(edge)>0) {
    if(edges=='omit') {
      fp <- fp[-edge,]
    } else {
      fp$edges[edge] <- T
    }
  }
  ## Detect (& possibly omit) outliers:
  if(nrow(fp)>2) {
    out <- outlier(fp[,1], logical=T)
    if(any(out)) {
      if(which(out)!=which.max(fp$power)) {
        if(outliers=='omit') {
          fp <- fp[-which(out),]
        } else {
          fp$outliers <- out
        }
      }
    }
  }

  ## Peaks in reconstructed signal:
  frp <- findpeaks(r.sig)

  frp.p <- lapply(fp$peak, function(x) tail(which((x-frp[,2])>=0), 1))
  frp.p <- unlist(lapply(frp.p, function(x) {
    if(length(x)==0) {
      NA
    } else {
      x
    }}))

  fp$rs.pk <- frp[frp.p,2]

  frp <- frp[frp.p,]
  if(class(frp)!='matrix') frp <- matrix(frp, ncol=4)
  which.true <- which(apply(fp[,c(5,6)], 1, function(x) all(!x)))

  if(plotting) {
    par(mfrow=c(4,1), mar=c(2,5,1,1))
    plot(-dat$depth, type='n', xaxs='i', axes=F, ylab='Depth (m)')
    axis(1, labels=F)
    axis(2)
    box()

    if('lunge' %in% names(dat)) {
      abline(v=which(dat$lunge), lty=2, xpd=NA, col='grey')
      points(which(dat$lunge), -dat$depth[which(dat$lunge)],
             cex=2, pch=21, bg=3, xpd=NA)
    }
    lines(-dat$depth)

    plot(dat$d.depth, type='o', col='slategrey', xaxs='i',
         axes=F, ylab='Depth change')
    lines(dat$w.mn, col=2)
    axis(1, labels=F)
    axis(2)
    box()

    plot(w.pow, type='l', xaxs='i', axes=F, ylab='Weighted power')
    points(fp[,2], fp[,1], pch=21, bg='grey')

    for(i in 1:nrow(fp)) {
      if(!fp$outliers[i]& !fp$edges[i]) {
        lines(c(fp$start[i]:fp$end[i]), w.pow[c(fp$start[i]:fp$end[i])], col=4, lwd=3)
      }
    }

    abline(v=fp$start[which.true], col=4, lty=2)
    abline(v=fp$end[which.true], col=4, lty=2)

    points(fp[which.true,2], fp[which.true,1],
           pch=21, bg=4, cex=2, xpd=NA)
    axis(1, labels=F)
    axis(2)
    box()


  ##  par(cex.lab=1, cex.axis=1)
  ##  wt.image(dw, plot.legend = F, col.contour=1, plot.ridge=F,
  ##           siglvl=w.sig.lev)
    plot(r.sig, type='l', xaxs='i', ylab='Reconstructed\nsignal')
    points(frp[,2], frp[,1], pch=21, bg='grey')
    points(frp[which.true,2], frp[which.true,1],
           pch=21, bg=2, cex=2, xpd=NA)
  }

  if(return.wavelet) {
    list(parameters=parameters, wavelet=dw, peaks=fp)
  } else {
    fp
  }
}


#' \code{wt.image.plot} Plot wavelet image for single dive
#'
#' @param dw.obj Object (output from \code{d.wave}) for which wavelet image should be plotted
#' @return Returns a colour image showing the power across time (horizontal) and wavelet periods (vertical)
#' @details
#' @family Lunge detector functions
#' @seealso \code{\link{read.lunges}} to read lunge event timestamps,
#'   \code{\link{read.depth}} to read depth time series,
#'   \code{\link{join.dstamp}} to merge lunge events and dive data,
#'   \code{\link{number.dive}} to add dive numeric ID column,
#'   \code{\link{d.wave}} 2D lunge detector for individual dive
#' @author Martin Biuw
#' @examples
#' d <- read.depth(whalename='mn16_Jan25b')
#' t <- read.lunges(whalename='mn16_Jan25b')
#' df <- join.dstamp(d, t)
#' df$dive <- number.dive(df$depth)
#' d103 <- d.wave(dive=103)
#' plot.wt.image(d103)
#' @export

wt.image.plot <- function(dw.obj) {
  par(mfrow=c(1,1), mar=c(3,4,1,1))
  wt.image(dw.obj$wavelet, plot.contour=F, plot.ridge=F, graphics.reset=T)

  plot.levs <- dw.obj$parameters[[3]]

  plot.levs <- match(plot.levs, dw.obj$wavelet$Period)
  plot.levs <- plot.levs + (c(0,-dw.obj$parameters[[4]], dw.obj$parameters[[4]]))
  plot.levs <- seq(par('usr')[3], par('usr')[4],
                   length=length(dw.obj$wavelet$Period))[plot.levs]
  segments(rep(par('usr')[1], length(plot.levs)), plot.levs,
           rep(par('usr')[2]*0.85, length(plot.levs)),
           lty=c(1,2,2))
  pks <- (dw.obj$peaks$rs.pk[which(!dw.obj$peaks$outlier & !dw.obj$peaks$edge)]/
    dim(dw.obj$wavelet$Power)[2])*0.85

  for(i in 1:length(pks)) {
    lines(rep(pks[i], 2), plot.levs[c(2,3)], lty=3)
  }
}

## #' \code{match.lunges} Deprecated lunge matching function
## #' @param
## #' @return
## #' @details
## #' @family
## #' @seealso
## #' @author Martin Biuw
## #' @examples
## #' @export

## match.lunges <- function(data=df) {
##   ## Add option for ignoring outliers and edge points in matching
##
##   true.pos <- data$lunge & data$pk.reg
##   false.neg <- data$lunge & !data$pk.reg
##
##   true.lunges <- which(data$lunge)
##   len.pk.reg <- round(mean(rle(data$pk.reg)$lengths[which(rle(data$pk.reg)$values)])/2)
##
##   lunge.reg <- unlist(lapply(true.lunges, function(x) {
##     c((x-len.pk.reg):(x+len.pk.reg))
##   }))
##
##   lunge.regs <- rep(F, nrow(data))
##   lunge.regs[lunge.reg] <- T
##
##   false.pos <- !lunge.regs & data$pks
##   true.neg <- !lunge.regs & !data$pk.reg
##
##   conf.mat <- matrix(c(length(which(true.pos)), length(which(false.neg)),
##                        length(which(false.pos)), length(which(true.neg))),
##                      ncol=2, byrow=T)
##   dimnames(conf.mat) <- list(Observed=c('True', 'False'), Predicted=c('True', 'False'))
##
##
## ##  tpr <- conf.mat[1,1]/apply(conf.mat, 2, sum)[1]
## ##  fnr <- conf.mat[2,1]/apply(conf.mat, 2, sum)[1]
## ##  fpr <- conf.mat[1,2]/apply(conf.mat, 2, sum)[2]
## ##  tnr <- conf.mat[2,2]/apply(conf.mat, 2, sum)[2]
##
##   list(matched=true.pos, confusion=conf.mat)
## }


