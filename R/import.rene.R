#' \code{matLunge} Read lunge data from acceleration data logger stored in Matlab format
#' @param file File name (including relative path) of file to be read.
#' @return Returns list with components representing lunges
#' @details
#' @family Lunge detector functions
#' @seealso \code{\link{matDive}} to read dive data,
#' \code{ncWhale} to read data from standard netCDF tagtools file,
#' \code{gpsWhale} to read GPS position data,
#' \code{whale2dstamp} to convert to dstamp format that can be used in 2D lunge detection with validation (e.g. \code{run.d.wave})
#' @author Martin Biuw
#' @examples
#' @export

matLunge <- function(file='./Data/Rene data/lunge_res/03_Herring_mn13_340a_lunge_herrspeed.mat') {
  require(R.matlab)
  mat <- readMat(file)$lunge
  names(mat) <- c("depid", "method", "size", "n", "ndives", "diveID",
                  "npeak", "start", "end", "tmax", "vmax", "percent",
                  "th", "tblank", "phase", "phid", "tdive", "true_time",
                  "true_max", "true_dive_id", "true_phase_txt", "true_phase_id",
                  "true_descent_speed", "ds")

  mat <- mat[-which(names(mat) %in%
                      c('phase', 'phid', 'true_max',
                        'true_phase_txt', 'true_phase_id'))]

  classes <- unlist(lapply(mat, class))
  lengths <- as.numeric(summary(mat)[,1])
  which.simplify <- which(lengths==1)
  for(i in which.simplify) {
    if(!is.null(mat[[i]])) {
      mat[[i]] <- as.vector(mat[[i]])
    }
  }

  dims <- lapply(mat, function(x) dim(unlist(x)))
  which.simplify <- which(unlist(lapply(dims, length))==2)
  for(i in which.simplify) {
    if(any(dim(mat[[i]])==1)) {
      mat[[i]] <- as.vector(mat[[i]])
    }
  }

  mat
}

#' \code{matDive} Read dive data from acceleration data logger stored in Matlab format
#' @param file File name (including relative path) of file to be read.
#' @return Returns list with components representing dives
#' @details
#' @family Lunge detector functions
#' @seealso \code{\link{matLunge}} to read lunge data,
#' \code{ncWhale} to read data from standard netCDF tagtools file,
#' \code{gpsWhale} to read GPS position data,
#' \code{whale2dstamp} to convert to dstamp format that can be used in 2D lunge detection with validation (e.g. \code{run.d.wave})
#' @author Martin Biuw
#' @examples
#' @export

matDive <- function(file='./Data/Rene data/lunge_res/03_Herring_mn13_340a_herrdive.mat') {
  mat <- readMat(file)$dive
  names(mat) <- c("depid", "method", "units", "sampling_rate", "start",
                  "end", "max", "cue_max", "duration", "descent_start",
                  "descent_end", "descent_duration", "descent_gradient",
                  "descent_vertvel_minmax", "descent_vertvel_mu",
                  "descent_vertvel_prctile", "descent_vertacc_minmax",
                  "descent_vertacc_mu", "descent_vertacc_prctile",
                  "descent_ocdr_minmax", "descent_ocdr_mu",
                  "descent_ocdr_prctile", "descent_paddle_minmax",
                  "descent_paddle_mu", "descent_paddle_prctile",
                  "descent_kalman_minmax", "descent_kalman_mu",
                  "descent_kalman_prctile", "descent_jerk_minmax",
                  "descent_jerk_mu", "descent_jerk_prctile", "bottom_start",
                  "bottom_end", "bottom_duration", "bottom_depth_start",
                  "bottom_depth_end", "bottom_depth_minmax", "bottom_depth_mu",
                  "bottom_depth_prctile", "bottom_depth_gradient",
                  "bottom_vertvel_minmax", "bottom_vertvel_mu",
                  "bottom_vertvel_prctile", "bottom_vertacc_minmax",
                  "bottom_vertacc_mu", "bottom_vertacc_prctile",
                  "bottom_ocdr_minmax", "bottom_ocdr_mu", "bottom_ocdr_prctile",
                  "bottom_paddle_minmax", "bottom_paddle_mu", "bottom_paddle_prctile",
                  "bottom_kalman_minmax", "bottom_kalman_mu", "bottom_kalman_prctile",
                  "bottom_jerk_minmax", "bottom_jerk_mu", "bottom_jerk_prctile",
                  "ascent_start", "ascent_end", "ascent_duration", "ascent_gradient",
                  "ascent_vertvel_minmax", "ascent_vertvel_mu", "ascent_vertvel_prctile",
                  "ascent_vertacc_minmax", "ascent_vertacc_mu", "ascent_vertacc_prctile",
                  "ascent_ocdr_minmax", "ascent_ocdr_mu", "ascent_ocdr_prctile",
                  "ascent_paddle_minmax", "ascent_paddle_mu", "ascent_paddle_prctile",
                  "ascent_kalman_minmax", "ascent_kalman_mu", "ascent_kalman_prctile",
                  "ascent_jerk_minmax", "ascent_jerk_mu", "ascent_jerk_prctile",
                  "surface_start", "surface_end", "surface_duration", "surface_blows",
                  "surface_interval_minmax", "surface_interval_mu")
  mat
}


#' \code{ncWhale} Read all data from acceleration data logger stored in standard netCDF format compatible with \code{tagtools}
#' @param file File name (including relative path) of file to be read.
#' @return Returns list with components representing dives
#' @details
#' @family Lunge detector functions
#' @seealso \code{\link{matLunge}} to read lunge data,
#' \code{matDive} to read dive data,
#' \code{gpsWhale} to read GPS position data,
#' \code{whale2dstamp} to convert to dstamp format that can be used in 2D lunge detection with validation (e.g. \code{run.d.wave})
#' @author Martin Biuw
#' @examples
#' @export

ncWhale <- function(file='./Data/Rene data/lunge_res/03_Herring_mn13_340a_Herr.nc') {
  require(ncdf4)
  nc <- nc_open(file)
  nc.atts <- ncatt_get(nc, 0)
  nc.names <- names(nc$var)
  out.list <- vector(mode='list')
  out.list$globals <- nc.atts
  for(i in 1:length(nc.names)) {
    out.list[[i+1]] <- ncvar_get(nc, nc.names[i])
  }
  nc_close(nc)
  names(out.list)[-1] <- nc.names
  out.list
}

#' \code{gpsWhale} Read Fastloc gps data from ESRI shapefile.
#' @param file File name (including relative path) of file to be read.
#' @param crs Coordinate Reference System (CRS) to be used (default is UTM zone 33, appropriate for Norwegian coast).
#' @return Returns SpatialPointsDataFrame
#' @details
#' @family Lunge detector functions
#' @seealso \code{\link{matLunge}} to read lunge data,
#' \code{\link{matDive}} to read dive data,
#' \code{ncWhale} to read data from standard netCDF tagtools file,
#' \code{whale2dstamp} to convert to dstamp format that can be used in 2D lunge detection with validation (e.g. \code{run.d.wave})
#' @author Martin Biuw
#' @examples
#' @export

gpsWhale <- function(filename='./Data/Rene data/shp/Mn-2014-Dec-16.shp',
                     crs='+proj=utm +zone=33 +ellps=intl +units=m +no_defs') {
  require(rgdal)
  gps <- readOGR(filename)
  gps$Time <- as.POSIXct(strptime(as.character(gps$Time),
                                  '%d/%m/%Y %H:%M:%S'), tz='GMT')
  if(!is.na(crs)) {
    gps <- spTransform(gps, CRS(crs))
  }
  gps
}


