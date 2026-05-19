\name{create.ctdata}
\alias{create.ctdata}
\alias{join.ctdata}
\title{Create and combine MONOCAR data sets}
\description{Create and combine data sets for estimation of MONOCAR
  models by \code{monocar.estimate} and \code{monocar.hist}}
\usage{
create.ctdata(x, v, t1, t2 = NULL,
              data = NULL,
              series.name = NULL,
              house.name = NULL,
              process.name = NULL,
              transform = NULL,
              transform.var = NULL,
              time.origin = NULL,
              time.unit = NULL,
              exo.data = NULL,
              inclusive.end.date = FALSE)

join.ctdata(\dots)
}
\arguments{
  \item{x}{A numeric vector giving the value of each observation.}
  \item{v}{A numeric vectorgiving the (untransformed) variance of each
    observation.}
  \item{t1}{A numeric, date-time, or date vector giving
    the start time for each observation.}
  \item{t2}{An optional numeric, date-time, or date vector
    giving the end time for each observation. Must be of the same class
    as \code{t1}. The default is to assume \code{t2} equals \code{t1},
    meaning instantaneous observations.}
  \item{data}{An optional data frame in which to search for variables.}
  \item{series.name}{An optional character, factor, or integer vector
    identifying the series for each observation. The default value of
    \code{NULL} uses the quoted \code{x} argument.}
  \item{house.name}{An optional character, factor, or integer vector
    identifying the house for each observation. The default is to
    use the series.}
  \item{process.name}{An optional character, factor, or integer vector
    identifying the process for each observation. The default is to
    assume all observations belong to the same process.}
  \item{transform}{An optional function specifying a transformation to
    apply to each observation. The function must take two inputs
    corresponding to vectors of means and variances. The default is
    leave the observations untransformed. Alternatively, a character
    object giving the name of a common variance-stabilizing
    transformation. The two implemented transformations are "arcsin" for
    binomial-distributed data (which uses
    \eqn{{arcsin}\left(\left(1-\frac{1}{5000}\right)\left(2\frac{x}{v}-1\right)\right)}{arcsin((1-1/5000)(2x/v-1))} and implied
    variance \eqn{\frac{1}{v}}{1/v}; this assumes \code{x} gives the
    number of successes and \code{v} gives the sample sizes for
    each observation) and "sqrt" for Poisson-distributed data (which
    uses \eqn{\sqrt{x+\frac{3}{8}}}{sqrt(x+3/8)} and implied variance
    \eqn{\frac{1}{4}}{1/4}; \code{v} is ignored).}
  \item{transform.var}{An optional function specifying a transformation
    to apply to the variance of each observation. The function must take
    two inputs corresponding to vectors of means and variances. The
    default is leave the variances untransformed or, if \code{transform}
    is a character object, the corresponding transformation of the
    variances for the named variance-stabilizing transformation.}
  \item{time.origin}{An optional numeric, date-time, or date object of
    length one specifying origin time to be subtracted from \code{t1}
    and \code{t2}.}
  \item{time.unit}{An optional numeric (if \code{t1} is a numeric
    vector) or difftime (if \code{t1} is a date-time or date vector)
    object specifying the units for \code{t1} and \code{t2} (i.e., the
    amount of time that should be considered to be length one.}
  \item{exo.data}{A vector, data frame, matrix, or formula of possible
    exogenous variables to be included in the data.}
  \item{inclusive.end.date}{If \code{t1} and \code{t2} are date vectors,
    whether the end date means the start or end of the day. Thus, for an
    observation which starts and ends on January 1, a \code{FALSE} value
    for \code{inclusive.end.date} would imply a length of zero hours
    while a \code{TRUE} value would imply a length of 24 hours.}
  \item{\dots}{Multiple arguments giving objects of class
  \code{ct.data.frame}, usually each created by \code{create.ctdata}, to
  be joined together. Alternatively, a single list containing multiple
  objects of class \code{ct.data.frame}.}
}
\details{
  See \code{\link{monocar.estimate}} for more details about the MONOCAR
  model.

  Observations are grouped into series, houses, and processes. Series
  correspond to different elements of the latent process, houses allow
  for house-specific effects, and processes allow for multiple,
  independent latent processes that share common parameters
  (specifically, \eqn{\Theta} and \eqn{\Sigma} as well as, optionally,
  \eqn{\mu} and/or \eqn{\delta}; whether the latter are included is
  controlled by \code{mu.by.process} and \code{delta.by.process}).
  Series are the most fundamental of the three. Houses are usually
  nested within series or identical to series. In most cases, all
  observations belong to the same process.

  If series are not specified, all observations will be assumed to
  belong to a single series. Usually this will not be the case. However,
  even when it is not the case, it may be easier to create separate
  \code{ct.data.frame} objects with \code{create.ctdata} for each series
  and combine them with \code{joint.ctdata}.
  
  If houses are not specified, they will be assumed to equal the series
  (implying series are not further divided into houses). If processes
  are not specified, all observations will be assumed to be based on a
  single (usually multivariate) latent process, as is typically the
  case.
}
\value{
  An object of class \code{ct.data.frame} containing the data.
}
\seealso{
  \code{\link{monocar.estimate}} and \code{\link{monocar.hist}}.
}
\examples{
## Create a simple data set with instantaneous observations every
## day for the first month of 2000
obs <- c(-0.25, 1.91, 0.06, 2.24, -0.48, -1.84, -1, -0.07,
-0.5, -1.17,  0.32, 0.67, 0.99, 0.39, 1.2, 0.12, -0.7, -0.65,
-1.52, -2.28,  0.98, 1.62, -0.36, 0.94, -1.61, -0.25, -0.76,
-1.81, -1.09, 0.9, -0.79)
obs.time <- seq(as.Date("2000-1-1"), by=1, length.out=length(obs))
obs.var <- rep(1, length(obs))
new.data <- create.ctdata(obs, obs.var, obs.time)
}
\keyword{monocar}
\keyword{ts}
