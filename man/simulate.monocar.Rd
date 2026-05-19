\name{simulate.monocar}
\alias{simulate.monocar}
\title{Simulate continuous-time data sets}
\description{Generate simulated continuous-time data sets}
\usage{
\method{simulate}{monocar}(object, nsim=1, seed=NULL,
  var=NULL, t1=seq(0,100,1.0), t2=NULL, data=NULL,
  series.name = NULL, house.name = NULL, process.name = NULL,
  transform = function(x, v) x, transform.var = function(x, v) v,
  exo.data=NULL, exovars=NULL, byexo=NULL, var.centers=NULL, \dots)
}
\arguments{
  \item{object}{A an object of class \code{monocar} or a list specifying
  the parameters of the model.}
  \item{nsim}{Number of data sets to simulate.  Defaults to \code{1}.}
  \item{seed}{An object specifying if and how the random number generator
    should be initialized. Either \code{NULL} or an integer that will be
    used in a call to \code{set.seed} before simulating the response
    vectors.  If set, the value is saved as the \code{seed}
    attribute of the returned value.  The default, \code{NULL} will
    not change the random generator state, and return
    \code{.Random.seed} as the \code{seed} attribute.}
  \item{var}{A numeric vector giving the variance of the measurement
    error for each observation. The default value of \code{NULL} means
    zero measurement error for each observation.}
  \item{t1}{A numeric vector giving the start time for each observation.}
  \item{t2}{An optional numeric vector giving the end time for each
    observation. By default, the end time of each observation equals the
    start time.}
  \item{data}{An optional data frame or object of class
  \code{ct.data.frame}.}
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
    apply to each observation. See \code{\link{create.ctdata}} for more
    details.}
  \item{transform.var}{An optional function specifying a transformation
    to apply to the variance of each observation. See
    \code{\link{monocar.data}} for more details.}
  \item{exo.data}{A vector, data frame, matrix, or formula of possible
    exogenous variables to be included in the data.}
  \item{exovars}{An optional formula giving the exogenous variables to
    be included, a character vector giving the names of exogenous
    variables to be included from \code{data}, or an integer vector
    giving the location of the exogenous variables to be included from
    \code{data}. The default is to include no exogenous variables.}
  \item{byexo}{A string or null value specifying whether exogenous
    variables are to have a single coefficient, separate coefficients
    for each series, or separate coefficients for each house.
    Possible values are \code{series}, \code{house} or
    \code{NULL}. Default is \code{NULL}, meaning a single coefficient,
    common to all series and houses, for each exogenous variable.}
  \item{var.centers}{For each house, the value around which to recenter
    the variances. See \code{\link{monocar.cal}} for more
    details.}
  \item{\dots}{Additional optional arguments.}
}
\value{
  If \code{nsim==1}, an object of class \code{ct.data.frame} containing
  the data. Otherwise, a list containing \code{nsim} objects of class
  \code{ct.data.frame}, each containing one simulated data set.
}
\seealso{
  \code{\link{create.ctdata}}, \code{\link{monocar.estimate}}, and \code{\link{monocar.hist}}.
}
\keyword{monocar}
\keyword{datagen}
\keyword{ts}
