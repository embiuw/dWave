#' \code{d.wave} 2D lunge detector for entire record, with validation
#' @param name Name of animal to analyse, if data are to be read from external files (default is NA)
#' @param data Data frame containing dive info (e.g. created using \code{read.depth}, \code{read.lunges},
#' \code{join.dstamp} and \code{number.dive}).
#' @param show.info Logical, should detaild information be shown while running function
#' @param plotting Whether to plot the results or not
#' @param match.window Size of time wondow within which peaks and lunges should be considered matched
#' @param smooth.window Initial window of smooth for gaussian running mean filter (equivalent to \code{presmooth.window} in \code{d.wave})
#' @param pper Peak wavelet period for weighted power analysis (equivalent to \code{pk.per} in \code{d.wave})
#' @param rper Range of all periods to be considered in weighted power analysis (equivalent to \code{per.rng} in \code{d.wave})
#' @param min.pdist Minimum period between consecutive peaks (equivalent to \code{min.pkdist} in \code{d.wave})
#' @param min.ppow Minimum amplitude for peak to be considered significant (equivalent to \code{a.pk} in \code{d.wave})
#' @param pk.e Proportion of a.pk, used to set the duration of what is considered a peak (equivalent to \code{pk.env} in \code{d.wave})
#' @return Returns a list with objects \code{parameters}, containing all parameters used when running detector,
#' \code{df} data frame containing original data plus peak detections and matching information, \code{detects} list containing wavelet ojects for each dive,
#' \code{matched} list containing 1) a logical vector indicating the occurrence of matched peaks/lunges, and 2) confusion matrix for matched and unmatched peaks/lunges
#' @details
#' @family Lunge detector functions
#' @seealso \code{\link{read.lunges}} to read lunge event timestamps,
#'   \code{\link{read.depth}} to read depth time series,
#'   \code{\link{join.dstamp}} to merge lunge events and dive data,
#'   \code{\link{number.dive}} to add dive numeric ID column,
#'   \code{\link{d.wave}} to run lunge detector omn single dives.
#' @author Martin Biuw
#' @examples
#' d <- read.depth(whalename='mn16_Jan25b')
#' t <- read.lunges(whalename='mn16_Jan25b')
#' df <- join.dstamp(d, t)
#' df$dive <- number.dive(df$depth)
#' df <- run.d.wave(data=df)
#' @export
#'

run.d.wave <- function(name=NA, data=NULL, show.info=F, plotting=T,
                       match.window=15, smooth.window=4, pper=16, rper=20,
                       min.pdist=40, min.ppow=0.05, pk.e=0.85) {

  used.args <- list(whalename=name, data=data,
                    presmooth.window=smooth.window,
                    pk.per=pper, per.rng=rper,
                    min.pkdist=min.pdist, a.pk=min.ppow, pk.env=pk.e,
                    return.wavelet=T, match.window=match.window)

  if(is.null(data)) {
    d <- read.depth(whalename=name)
    t <- read.lunges(whalename=name)
    df <- join.dstamp(d,t)
    df$dive <- number.dive(df$depth)
  } else {
    df <- data
  }

  detects <- list()
  df$edge <- df$out <- df$pk.power <- df$pks <- df$pk.reg <- rep(F, nrow(df))

  pb <- winProgressBar(title = "Running d.wave lunge detector:   ", min = 0,
                       max = max(df$dive, na.rm=T), width = 300)

  for(d in 1:max(df$dive, na.rm=T)) {
    setWinProgressBar(pb, d, title=paste("Running d.wave lunge detector:   ",
                                         round(d/max(df$dive, na.rm=T)*100, 0),
                                         "% done"))
    if(show.info) {
      cat('Running d.wave lunge detector for dive', d, 'of', max(df$dive, na.rm=T), '     \r')
      flush.console()
    }
    tmp <- try(d.wave(df, dive=d, presmooth.window=smooth.window,
                      pk.per=pper, per.rng=rper,
                      min.pkdist=min.pdist, a.pk=min.ppow, pk.env=pk.e,
                      return.wavelet=T, plotting=F),
               silent=T)
    detects[[d]] <- tmp
    if(class(tmp)!='try-error') {
      if(class(tmp$peaks)=='data.frame') {
        if(nrow(tmp$peaks)>=1) {
          df$pks[which(df$dive==d)[tmp$peaks$rs.pk]] <- T
          df$pk.power[which(df$dive==d)[tmp$peaks$peak]] <- tmp$peaks$power
          df$out[which(df$dive==d)[tmp$peaks$peak]] <- tmp$peaks$outliers
          df$edge[which(df$dive==d)[tmp$peaks$peak]] <- tmp$peaks$edges
          for(p in 1:nrow(tmp$peaks)) {
            pr <- c(tmp$peaks$start[p]:tmp$peaks$end[p])
            df$pk.reg[which(df$dive==d)[pr]] <- T
          }
        }
      }
    }
  }

  close(pb)

  if(show.info) {
    cat('\n\nDone!\n\n')
    flush.console()
  }

  if('lunge' %in% names(df)) {
    matched <- match.lunges(data=df)
    df$matched <- rep(FALSE, nrow(df))
    df$matched[which(df$pks)] <- matched$matched
  }

  if(plotting) {
    par(mfrow=c(1,1))
    plot(-depth~time, data=df, type='l', axes=F)
    axis.POSIXct(1, at=pretty(df$time, length=5), format='%H:%M')
    axis(2)
    text(as.numeric(aggregate(df$time, list(df$dive), mean)$x),
         rep(par('usr')[4], max(df$dive, na.rm=T)),
         c(1:max(df$dive, na.rm=T)), pos=3, xpd=NA)
    if('lunge' %in% names(df)) {
      points(-depth~time, data=df[which(df$lunge),], col=2, cex=2)
    }
    points(-depth~time, data=df[which(df$pks),], col=3, pch=21, cex=2)
    if('lunge' %in% names(df)) {
      points(-depth~time, data=df[which(df$matched),], bg=rgb(0,1,0,0.5), col=2, pch=21, cex=2)
      print(matched$confusion)
    }
  }
  if('lunge' %in% names(df)) {
    list(parameters=used.args, df=df, detects=detects, matched=matched)
  } else {
    list(parameters=used.args, df=df, detects=detects)
  }
}


#' \code{optim.d.wave} Obtain ptimize
#' @param par Vector of starting values for parameters to be estimated
#' @param theData Data frame containing dive info (e.g. created using \code{read.depth}, \code{read.lunges},
#' \code{join.dstamp} and \code{number.dive}).
#' @return Returns the objective function for optimal parameter combinations
#' @details NOT YET FULLY IMPLEMENTED. EXTREMELY TIME CONSUMING!
#' @family Lunge detector functions
#' @seealso \code{\link{read.lunges}} to read lunge event timestamps,
#'   \code{\link{read.depth}} to read depth time series,
#'   \code{\link{join.dstamp}} to merge lunge events and dive data,
#'   \code{\link{number.dive}} to add dive numeric ID column,
#'   \code{\link{d.wave}} to run lunge detector omn single dives.
#'   \code{\link{run.d.wave}} to run \code{d.wave} on entire data record from an animal.
#' @author Martin Biuw
#' @examples
#' @export

optim.d.wave <- function(par, theData=df) {
  par[2] <- round(par[2])
  par[3] <- round(par[3])
  par[4] <- round(par[4])
  par[5] <- round(par[5])
  tmp <- try(run.d.wave(NA, theData, show.info=F, plotting=F, match.window=par[1], smooth.window=par[2], pper=par[3],
                    rper=par[4], min.pdist=par[5], min.ppow=par[6], pk.e=par[7]), silent=T)
  if(class(tmp) != 'try-error') {
    tpr <- tmp$matched$confusion[1,1]/apply(tmp$matched$confusion, 1, sum)[1]
    fnr <- tmp$matched$confusion[1,2]/apply(tmp$matched$confusion, 1, sum)[1]
    fpr <- tmp$matched$confusion[2,1]/apply(tmp$matched$confusion, 1, sum)[2]
    tnr <- tmp$matched$confusion[2,2]/apply(tmp$matched$confusion, 1, sum)[2]

    precision <- tmp$matched$confusion[1,1]/apply(tmp$matched$confusion, 2, sum)[1]

    lrp <- tpr/fpr
    lrn <- fnr/tnr
    dor <- lrp/lrn
    f.score <- 2 / ((1 / tpr) + (1 / precision))

    list(tpr=tpr, fnr=fnr, fpr=fpr, tnr=tnr, lrp=lrp, lrn=lrn, f.score=f.score)

    ##sum(as.vector(tmp$matched$confusion)[c(2,3)])/tmp$matched$confusion[1,1]
  } else {
    99999
  }
}


## ppers <- seq(4, 32, by=2)
## sensitivity <- specificity <- rep(NA, length(ppers))
## for(i in 1:length(ppers)) {
##   rwav <- optim.d.wave(par=c(4,ppers[i], 20, 40, 0.05, 0.85), theData=df)
##   sensitivity[i] <- rwav[1]
##   specificity[i] <- rwav[4]
## }

## Example:
## opt.tst <- nlminb(unlist(formals(run.d.wave)[-c(1:5)]),
##                   optim.d.wave,
##                   lower=c(2, 8, 4, 16, 20, 0, 0.5),
##                   upper=c(10, 32, 16, 128, 120, 0.2, 0.95))

## mn16_Jan19a.opt <- nlminb(unlist(formals(run.d.wave)[-c(1:4)]),
##                           optim.d.wave,
##                           upper=c(10, 32, 16, 128, 120, 0.2, 30),
##                           lower=c(2, 8, 4, 16, 20, 0, 5),
##                           theData=mn16_Jan19a,
##                           control=list(abs.tol=1e-20))

## mn16_Jan25a.opt <- nlminb(unlist(formals(run.d.wave)[-c(1:4)]),
##                           optim.d.wave,
##                           upper=c(10, 32, 16, 128, 120, 0.2, 30),
##                           lower=c(2, 8, 4, 16, 20, 0, 5),
##                           theData=mn16_Jan25a,
##                           control=list(abs.tol=1e-20))

## mn16_Jan25b.opt <- nlminb(unlist(formals(run.d.wave)[-c(1:4)]),
##                           optim.d.wave,
##                           upper=c(10, 32, 16, 128, 120, 0.2, 30),
##                           lower=c(2, 8, 4, 16, 20, 0, 5),
##                           theData=mn16_Jan25b,
##                           control=list(abs.tol=1e-20))

## mn16_Jan26.opt <- nlminb(unlist(formals(run.d.wave)[-c(1:4)]),
##                           optim.d.wave,
##                           upper=c(10, 32, 16, 128, 120, 0.2, 30),
##                           lower=c(2, 8, 4, 16, 20, 0, 5),
##                           theData=mn16_Jan26,
##                          control=list(abs.tol=1e-20))

## mn17_022LLa.opt <- nlminb(unlist(formals(run.d.wave)[-c(1:4)]),
##                          optim.d.wave,
##                          upper=c(10, 32, 16, 128, 120, 0.2, 30),
##                          lower=c(2, 8, 4, 16, 20, 0, 5),
##                          theData=mn17_022LLa,
##                          control=list(abs.tol=1e-20))

## mn17_022LLb.opt <- nlminb(unlist(formals(run.d.wave)[-c(1:4)]),
##                           optim.d.wave,
##                           upper=c(10, 32, 16, 128, 120, 0.2, 30),
##                           lower=c(2, 8, 4, 16, 20, 0, 5),
##                           theData=mn17_022LLb,
##                           control=list(abs.tol=1e-20))

## mn17_026LLa.opt <- nlminb(unlist(formals(run.d.wave)[-c(1:4)]),
##                           optim.d.wave,
##                           upper=c(10, 32, 16, 128, 120, 0.2, 30),
##                           lower=c(2, 8, 4, 16, 20, 0, 5),
##                           theData=mn17_026LLa,
##                           control=list(abs.tol=1e-20))

## mn18_013LLa.opt <- nlminb(unlist(formals(run.d.wave)[-c(1:4)]),
##                           optim.d.wave,
##                           upper=c(10, 32, 16, 128, 120, 0.2, 30),
##                           lower=c(2, 8, 4, 16, 20, 0, 5),
##                           theData=mn18_013LLa,
##                           control=list(abs.tol=1e-20))

## all.pars <- as.data.frame(rbind(mn16_Jan19a.opt$par,
## mn16_Jan25a.opt$par,
## mn16_Jan25b.opt$par,
## mn16_Jan26.opt$par,
## mn17_022LLa.opt$par,
## mn17_022LLb.opt$par,
## mn17_026LLa.opt$par,
## mn18_013LLa.opt$par))
## row.names(all.pars) <- whalenames

