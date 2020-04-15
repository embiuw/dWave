#' \code{plot.TDR.detects} Plot results from runing \code{d.wave.TDR}
#' @param data Data list containing output from \code{d.wave.TDR}.
#' @param pType Point type for lunge detections
#' @param pCol Point colour for lunge detections
#' @param pSize Point size for lunge detections
#' @param pAlpha Point transparency for lunge detections
#' @param Title Title of pot (usually used to give name of the animal)
#' @family Lunge detector functions
#' @seealso \code{\link{read.lunges}} to read lunge event timestamps,
#'   \code{\link{read.depth}} to read depth time series,
#'   \code{\link{join.dstamp}} to merge lunge events and dive data,
#'   \code{\link{number.dive}} to add dive numeric ID column,
#'   \code{\link{d.wave}} to run lunge detector omn single dives.
#'   \code{\link{d.wave.TDR}} to run lunge detector on entire TDR record.
#' @author Martin Biuw
#' @examples
#' tst <- d.wave.TDR()
#' plot.TDR.detects(tst)
#' par(mfrow=c(2,2), mar=c(2,4,1,1))
#' plot.TDR.detects(tst)
#' plot.TDR.detects(pType='+', pCol=2, pSize=1)
#' plot.TDR.detects(pType='+', pCol=4, pSize=1, pAlpha=1)
#' plot.TDR.detects(pType=19, pCol=3, pAlpha=1)
#' @export
#'

plot.TDR.detects <- function(data=tst, pType=21, pCol=3, pSize=2, pAlpha=0.5, Title='Test') {
  df <- tst$df
  plot(depth~time, data=df, type='l', axes=F, col='grey',
       ylim=c(max(df$depth, na.rm=T), 0))
  axis.POSIXct(1, at=pretty(df$time, length=5), format='%H:%M')
  axis(2)
  abline(v=par('usr')[1])
  abline(h=par('usr')[3])
  abline(h=0, col='grey')
  if(!is.na(Title)) text(mean(par('usr')[c(1,2)]), par('usr')[4], Title, adj=c(0.5,0), xpd=NA)

  pCol <- as.vector(col2rgb(pCol))/255
  pCol <- rgb(pCol[1], pCol[2], pCol[3], pAlpha)

  if('pks' %in% names(df)) {
    if(pType %in% c(21:25)) {
      points(depth~time, data=df[which(df$pks),], pch=pType, bg=pCol, cex=pSize)
    } else {
      points(depth~time, data=df[which(df$pks),], pch=pType, col=pCol, cex=pSize)
    }
  }
}


