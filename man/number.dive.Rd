% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detect.lunges.R
\name{number.dive}
\alias{number.dive}
\title{\code{number.dive} Split data into distinct dives}
\usage{
number.dive(data = df$depth, d.cut = 5)
}
\arguments{
\item{data}{Depth variable in data frame created using \code{join.dstamp}}

\item{d.cut}{Initial cutoff (in meters) for dives to be detected.}
}
\value{
Returns a data frame with variables \code{time},
\code{depth} and \code{lunge}
}
\description{
\code{number.dive} Split data into distinct dives
}
\examples{
t <- number.dive(df)
}
\seealso{
\code{\link{read.lunges}} to read lunge event timestamps,
  \code{\link{read.depth}} to read dive data,
  \code{\link{join.dstamp}} to add dive numeric ID column,

Other Lunge detector functions: \code{\link{d.wave.TDR.original}},
  \code{\link{d.wave.TDR}}, \code{\link{d.wave}},
  \code{\link{dstamp.plot}}, \code{\link{filter.lunges}},
  \code{\link{gpsWhale}}, \code{\link{import.all.rene}},
  \code{\link{join.dstamp}}, \code{\link{matDive}},
  \code{\link{matLunge}}, \code{\link{match.lunges}},
  \code{\link{mkDygraph}}, \code{\link{ncWhale}},
  \code{\link{optim.d.wave}}, \code{\link{read.depth}},
  \code{\link{read.lunges}}, \code{\link{run.d.wave}},
  \code{\link{whale2dstamp}}, \code{\link{wt.image.plot}}
}
\author{
Martin Biuw
}
\concept{Lunge detector functions}
