#' \code{import.all.rene} Import all data records from Dtag /Little Leonardo deployments
#' @param directory Directory where dive /accelerometer data are stored
#' @param gps.dir Directory where GPS data are stored
#' @param crs Coordinate Reference System (CRS) to be used (default is UTM zone 33, appropriate for Norwegian coast).
#' @return Returns a list with components \code{raw}, containing the raw imported dive/accelerometer data,
#' \code{dive} containing specific dive data, \code{lunge} containing lunge data
#' @details NOTE! This function creates new objects in the working environment automatically,
#' with names corresponding to the whale deployment names. No need therefore to save the results
#' into a new object.
#' @family Lunge detector functions
#' @seealso \code{\link{matLunge}} to read lunge data,
#' \code{\link{matDive}} to read dive data,
#' \code{ncWhale} to read data from standard netCDF tagtools file,
#' \code{gpsWhale} to read GPS position data,
#' \code{whale2dstamp} to convert to dstamp format that can be used in 2D lunge detection with validation (e.g. \code{run.d.wave})
#' @author Martin Biuw
#' @examples
#' import.all.rene()
#' ## Then add new object \code{ds} to each resulting list, in dstamp format ready for d.wave:
#' file.pos <- ls()[grep('mn1', ls())]
#' for(i in file.pos) eval(parse(text=paste(i, '$ds <- whale2dstamp(', i, ')', sep='')))
#' @export

import.all.rene <- function(directory='./Data/Rene data/lunge_res',
                            gps.dir='./Data/Rene data/shp',
                            gpsCRS='+proj=utm +zone=33 +ellps=intl +units=m +no_defs') {
  whalefiles <- dir(directory)
  wframes <- whalefiles[grep('.nc', whalefiles)]
  lungefiles <- whalefiles[grep('lunge', whalefiles)]
  divefiles <- whalefiles[grep('dive', whalefiles)]

  whalenames <- unlist(lapply(wframes, function(x) {
    gsub('.nc', '', paste(unlist(strsplit(x, '_'))[3],
                          unlist(strsplit(x, '_'))[4], sep='_'))
  }))
  whalenames <- c(whalenames, unlist(lapply(lungefiles, function(x) {
    gsub('.mat', '', paste(unlist(strsplit(x, '_'))[3],
                           unlist(strsplit(x, '_'))[4], sep='_'))
  })))
  whalenames <- c(whalenames, unlist(lapply(divefiles, function(x) {
    gsub('.mat', '', paste(unlist(strsplit(x, '_'))[3],
                           unlist(strsplit(x, '_'))[4], sep='_'))
  })))

  whalenames <- unique(whalenames)

  for(i in 1:length(whalenames)) {
    cat('Reading data for', whalenames[i], '(', i, 'of', length(whalenames), ')\n')
    flush.console()
    r <- grep(whalenames[i], wframes)
    raw.str <- paste('ncWhale("', paste(directory, wframes[r], sep='/'), '")', sep='')

    d <- grep(whalenames[i], divefiles)
    dive.str <- paste('matDive("', paste(directory, divefiles[d], sep='/'), '")', sep='')

    l <- grep(whalenames[i], lungefiles)
    lunge.str <- paste('matLunge("', paste(directory, lungefiles[l], sep='/'), '")', sep='')

    lst <- eval(parse(text=paste('lst <- list(raw=', raw.str, ', dive=',
                                 dive.str, ', lunge=', lunge.str, ')')))
    if(!is.na(gps.dir)) {
      dat <- lst$raw$globals$dephist_deploy_datetime_start
      ddat <- as.POSIXct(strptime(dat, '%d-%b-%Y %H:%M:%S'), tz='GMT')
      if(is.na(ddat)) {
        ddat <- as.POSIXct(strptime(dat, '%Y-%m-%d %H:%M:%S'), tz='GMT')
      }
      fdat <- format(ddat, '%Y-%b-%d')
      which.gps <- grep(paste(fdat, 'shp', sep='.'), dir(gps.dir))
      if(length(which.gps)==1) {
        lst$gps <- gpsWhale(paste(gps.dir, dir(gps.dir)[which.gps], sep='/'))
        if(proj4string(lst$gps)!=gpsCRS) {
          lst$gps <- spTransform(lst$gps, CRS(gpsCRS))
        }
      }
    }
    eval(parse(text=paste('assign("', whalenames[i],
                          '", lst, .GlobalEnv)',sep='')))
  }
}


#' \code{whale2dstamp} Convert data from object created by \code{import.all.rene} to data frame compatible with \code{d.wave} and \code{run.d.wave}.
#' @param lst Object obtained from \code{ncWhale}
#' @return Returns data frame that can be used in 2D lunge detector with validation
#' @details
#' @family Lunge detector functions
#' @seealso \code{import.all.rene},
#' \code{\link{matLunge}} to read lunge data,
#' \code{\link{matDive}} to read dive data,
#' \code{ncWhale} to read data from standard netCDF tagtools file,
#' \code{gpsWhale} To read Fastloc GPS data
#' @author Martin Biuw
#' @examples
#' @export

whale2dstamp <- function(lst=mn13_340a) {
  d <- lst$raw$P
  t <- lst$raw$globals$dephist_deploy_datetime_start
  tt <- as.POSIXct(strptime(t, '%d-%b-%Y %H:%M:%S'), tz='UTC')
  if(is.na(tt)) {
    tt <- as.POSIXct(strptime(t, '%Y-%m-%d %H:%M:%S'), tz='UTC')
  }
  t <- tt+(c(0:(length(d)-1)))
  lunge <- rep(F, length(d))
  lunge.t <- lst$lunge$true_time
  lunge.t <- lunge.t[which(!is.nan(lunge.t))]
  lunge[lunge.t] <- T
  dive <- rep(NA, length(d))
  dive.per <- cbind(lst$dive$start, lst$dive$end)
  for(i in 1:nrow(dive.per)) {
    dive[c(dive.per[i,1]:dive.per[i,2])] <- i

  }
  data.frame(time=t, depth=d, lunge=lunge, dive=dive)
}


