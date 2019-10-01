#' \code{d.wave.TDR} 2D lunge detector for individual dive in a TDR record
#' @param tdr.file R object containing pre-processed dive record. If set to NA, this function can run the detector on data from all individuals with data records stored in the directory specified in \code{tdr.folder}.
#' @param tdr.folder Folder in which TDR files are stored
#' @details Run the 2D lunged detector on TDR-type dive data, i.e. this version has no valiudation against pre-determined lunges
#' @family Lunge detector functions
#' @seealso \code{\link{run.d.wave}} for data type with validation against pre-determined lunges.
#' @author Martin Biuw
#' @examples
#' tdr.list <- d.wave.TDR(tdr.file=NA)
#' @export

d.wave.TDR <- function(tdr.file='./Data/HVTag/Rfiles/Whale_2015Feb23.RData',
                       tdr.folder='./Data/HVTag/Rfiles') {
  if(!is.na(tdr.file)) {
    load(tdr.file)
    whalename <- tail(unlist(strsplit(tdr.file, '/')),1)
    whalename <- gsub('.RData', '', whalename)
    eval(parse(text=paste('tdr <- ', whalename, '$original')))
    cat('Running 2D lunge detector for whale', whalename, '     \n')
    flush.console()
    tdr <- run.d.wave(data=tdr)
    tdr
  } else {
    tdr.file <- paste(tdr.folder, dir(tdr.folder), sep='/')
    tdr.list <- vector('list', length(dir(tdr.folder)))
    names(tdr.list) <- gsub('.RData', '', dir(tdr.folder))
    ## Note! Names ended up inconsistent from w=32 onwards.
    ## Fixed it externally, but check what's wrong here.

    for(w in 1:length(tdr.file)) {
      load(tdr.file[w])
      whalename <- tail(unlist(strsplit(tdr.file[w], '/')),1)
      whalename <- gsub('.RData', '', whalename)
      eval(parse(text=paste('tdr <- ', whalename)))
      cat('Running 2D lunge detector for whale',
          whalename, '(file', w, 'of', length(tdr.file), ')     \r')
      flush.console()
      tdr <- try(run.d.wave(data=tdr$original), silent=T)
      if(class(tdr)!='try-error') {
##        tdr.tmp$original <- tdr$df
##        dg <- mkDygraph(tdr)
##        tdr.tmp$graph <- dg
##        tdr.list[[w]] <- tdr.tmp
          tdr.list[[w]] <- tdr
        ##eval(parse(text=paste(whalename, '<- tdr.tmp')))
      }
    }
    tdr.list
  }
}

#' \code{filter.lunges} Filter out lunges during short & shallow dives
#' @param data Output from running \code{d.wave.TDR}.
#' @param duration Maximum dive duration of dives considered for filtering out lunges
#' @param depth Maximum dive depth of dives considered for filtering out lunges
#' @details
#' @family Lunge detector functions
#' @seealso \code{\link{d.wave.TDR}} to run 2D detector on TDR data.
#' @author Martin Biuw
#' @examples
#' tdr.list <- d.wave.TDR(tdr.file=NA)
#' @export

filter.lunges <- function(data=tdr.df, duration=60, depth=10) {
  data$bigdive <- rep(T, nrow(data))
  data$whale.dive <- paste(data$WhaleID, data$dive, sep='.')

  data$whale.dive[which(is.na(data$dive))] <- NA

  dive.depths <- aggregate(data$depth, list(data$whale.dive), max)
  dive.durs <- aggregate(data$depth, list(data$whale.dive), length)

  which.flag <- intersect(which(dive.depths$x<depth),
                          which(dive.durs$x<duration))

  pb <- winProgressBar(title = "Running lunge filter:   ", min = 0,
                       max = length(which.flag), width = 300)

  for(i in which.flag) {
    setWinProgressBar(pb, i, title=paste("Running lunge filter:   ",
                                         round(i/length(which.flag)*100, 0),
                                         "% done"))
    data$bigdive[which(data$whale.dive==dive.depths$Group.1[i])] <- F
  }
  close(pb)

  data[,-length(data)]
}
