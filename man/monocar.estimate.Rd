\name{monocar.estimate}
\alias{monocar.estimate}
\title{Estimate MONOCAR model}
\description{
  Estimate a Multivariate overlapping noisy observation continuous-time
  autoregressive (MONOCAR) model
}
\usage{
monocar.estimate(data, init=list(), restrict=list(),
    subset=NULL, exovars=NULL, byexo="series",
    enforce.bounds=TRUE, lower.bounds=NULL, upper.bounds=NULL,
    delta.group=TRUE, verbose=1, mu.by.process=FALSE,
    delta.by.process=FALSE, var.by.house=FALSE, remove.unused=TRUE,
    \dots)
}
\arguments{
  \item{data}{A data file of class \code{ct.data.frame} (see
    \code{\link{create.ctdata}})}
  \item{init}{A list containing initial values for parameters.
    Omitted items will be guessed. The default is for \code{theta},
    \code{sigma}, and \code{alpha} to be identity
    matrices, one of \code{mu} or \code{delta} to be the mean of each
    series or house (depending on whether \code{delta.group} is used),
    and variance-transformation parameters to be
    \code{0} for \code{var.add} and \code{1} for
    \code{var.mult} and \code{var.pow} (i.e., no transformation).}
  \item{restrict}{A list of which parameters to be restricted to
    their initial values. \code{TRUE} means restricted to initial
    value, \code{FALSE} means unrestricted, \code{NA} means
    restricted to equal the previous parameter (useful, for example, to
    set the variance-transformation parameters to be equal for multiple
    series. The parameter \code{sigma} can also be specified as one of
    the following strings: "diagonal" (meaning only the off-diagonal
    elements of \code{sigma} are restricted), "restricted" (meaning all
    elements of \code{sigma} are restricted), or "unrestricted" (meaning
    all elements of \code{sigma} are unrestricted). The default is for
    \code{theta} to be unrestricted, \code{sigma} to be restricted to be
    diagonal, \code{alpha} to be fully restricted, one of \code{mu} or
    \code{delta} to be unrestricted (depending on \code{delta.group} or
    whether initial values are explicitly specified), and
    variance-transformation parameters to be restricted unless initial
    values are explicitly specified.}
  \item{subset}{An optional vector specifying a subset of observations to be
    used in the fitting process.}
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
  \item{enforce.bounds}{A boolean or integer specifying whether to
    enforce bounds on the variance-transformation parameters during
    optimization. If equal to \code{FALSE} or \code{0}, no bounds are
    imposed. If equal to \code{TRUE} or \code{1}, \code{var.add} and
    \code{var.mult} are restricted to be non-negative, avoiding negative
    transformed variances. If equal to \code{2}, \code{var.pow} is also
    restricted to be non-negative, ensuring the variance transformation
    is non-decreasing. The default is \code{TRUE}.}
  \item{lower.bounds}{An optional vector of lower bounds for parameters if
    desired. The default is NULL, meaning no lower bounds other than
    those implied by \code{enforce.bounds}.}
  \item{upper.bounds}{An optional vector of lower bounds for parameters if
    desired. The default is NULL, meaning no lower bounds.}
  \item{delta.group}{An optional vector of integers specifying which
    houses are associated with which series. Alternatively, the logical
    value \code{TRUE}, in which case the data will be used to determine
    which series corresponds to each house, or
    \code{NULL}, in which case no attempt to recenter \eqn{\delta} will
    be taken. \code{delta.group} is used to recenter
    \eqn{\delta} such that the mean of \eqn{delta} within each series is
    zero with \eqn{\mu} adjusted to compensate (thus providing an
    equivalent parameterization). This only makes sense if houses are
    nested within series. The value \code{TRUE} will be treated as
    \code{NULL} if houses are not
    nested within series. The default is \code{TRUE}.}
  \item{verbose}{An integer specifying the amount of output. \code{0}
    means no output, \code{1} means a dot (\code{.}) for each gradient
    computation and a star (\code{*}) for each Hessian evaluation, and a
    \code{2} means a dot for each likelihood evaluation, a colon
    (\code{:}) for each gradient evaluation, and a star for each Hessian
    evaluation. A \code{+} indicates an optimization algorithm converged
    while \code{-} means it failed to do so. The default is \code{1}.}
  \item{mu.by.process}{If multiple processes are specified, a boolean
    indicating whether the series means should be identical across
    processes. The default is \code{FALSE}.}
  \item{delta.by.process}{If multiple processes are specified, a boolean
    indicating whether the house means should be identical across
    processes. The default is \code{FALSE}.}
  \item{var.by.house}{Whether variance-transformation parameters should
    be applied by house or by series. The default is to apply by
    series.}
  \item{remove.unused}{Whether to remove unused series and houses,
    thereby reducing the number of parameters. This can be particularly
    useful when using a subset of the data that does not contain all
    houses. Note that the size of the parameters specified in
    \code{init} and \code{restrict} may change as a result.}
  \item{\dots}{Additional parameters to be passed to the
    \code{monocar.cal}, which carries out the optimization.}
}
\details{
  Observations are grouped into series, houses, and processes. Series
  correspond to different elements of the latent process, houses allow
  for house-specific effects, and processes allow for multiple,
  independent latent processes that share common parameters
  (specifically, \eqn{\Theta} and \eqn{\Sigma} as well as, optionally,
  \eqn{\mu} and/or \eqn{\delta}; whether the latter are included is
  controlled by \code{mu.by.process} and \code{delta.by.process}).
  Series are the most fundamental of the three. Houses are usually
  nested within series or identical to series. In most cases, all
  observations belong to a single process.
  
  The latent process, \eqn{x_t} is an Ornstein-Uhlenbeck process
  defined by the stochastic differential equation
  \deqn{x_t = \Theta\left(\mu - x_t\right)dt + \Sigma^{\frac{1}{2}} dW_t}{x_t = \Theta(\mu - x_t) dt + \Sigma^(1/2) dW_t}
  where \eqn{\Theta} describes how each element of the latent process
  affects itself and other elements of the latent process,
  \eqn{\Sigma} is the instantaneous variance-covariance matrix for the
  error process, \eqn{\mu} is the long-run mean of \eqn{\mu},
  and \eqn{W_t} is a multivariate Wiener process.

  The observation equations are defined by
  \deqn{y_i = \int_{r_i}^{s_i} x_{s_i,t} dt + \delta_{h_i} + \epsilon_i}
  when \eqn{q_i < r_i} or
  \deqn{y_i = x_{s_i,q_i} + \delta_{h_i} + \epsilon_i}
  when \eqn{q_i = r_i} with
  \deqn{\epsilon_i \sim \mathcal{N}\left(0, \gamma_{s_i} + \nu_{s_i}
  \left(\frac{v_i}{\bar{v}_{s_i}}\right)^{\xi_{s_i}}\right)}{\epsilon_i ~
  N(0, \gamma_{s_i} + \nu_{s_i} (v_i/\bar{v}_{s_i})^\xi_{s_i})}
  and independent across observations. Here, \eqn{y_i} is the value of
  observation \eqn{i}, \eqn{v_i} is the untransformed (or naive)
  variance of the error term for the observation, \eqn{\delta_j} is a
  house-specific effect for
  house \eqn{j}, \eqn{h_i} is the house of equation \eqn{i}, \eqn{q_i}
  is the start time for observation \eqn{i}, \eqn{r_i} is the end time
  for observation \eqn{i}, \eqn{\bar{v}_j} is the variance center for
  house \eqn{j}, and \eqn{s_i} and \eqn{h_i}
  index the series and house of observation \eqn{i}.
  \eqn{y_i}, \eqn{v_i}, \eqn{q_i}, \eqn{r_i}, \eqn{s_i}, and \eqn{h_i}
  are all known
  (entered as \code{x}, \code{v}, \code{t1}, and \code{t2}
  in \code{\link{create.ctdata}}).
  If \code{var.by.house} is \code{TRUE}, then \eqn{\gamma}, \eqn{\nu},
  \eqn{\xi}, and \eqn{\bar{v}} are subscripted by \eqn{h_i} rather than
  \eqn{s_i}.
  
  The variance centers (\eqn{\bar{v}}) are used to improve parameter
  estimation and interpretation. By default, they are the mean of the
  untransformed variances for observations in each house. These default
  to the mean variance of each series or house, but other values can be
  specified (see \code{\link{monocar.cal}}). Aside from possible
  computational issues, a change of variance centers does not affect
  the model once \eqn{\nu} and \eqn{\xi} are adjusted appropriately.
  
  For multiple processes, the latent process is replaced
  by multiple latex processes with common parameters \eqn{\Theta} and
  \eqn{\Sigma} defined by
  \deqn{x_{k,t} = \Theta\left(\mu_k - x_{k,t}\right)dt +
    \Sigma^{\frac{1}{2}}
    dW_{k,t}}{x_{k,t} = \Theta (\mu_k - x_{k,t}) dt + \Sigma^(1/2) dW_{k,t}}
  for all \eqn{k}, where \eqn{k} indexes the latent process. The
  observation equations also likewise adjusted by replacing \eqn{x_t} by
  \eqn{x_{k_i,t}}, where \eqn{k_i} is the latent process observed by
  observation eqn{i}.

  For a structural MONOCAR model, the latent process becomes
  \deqn{A x_t = \Theta\left(\mu - x_t\right)dt + \Sigma^{\frac{1}{2}}
    dW_t}{A x_t = \Theta (\mu - x_t)dt + \Sigma^(1/2) dW_t}
  where \eqn{A} is
  a matrix describing instantaneous, non-gradual effects and is assumed
  to be the identity matrix by default.

  In specifying \code{init} and \code{restrict}, the possible parameters
  to list are
  \code{theta} (\eqn{\Theta}, an \eqn{N \times N} matrix), \code{sigma}
  (\eqn{\Sigma}, an \eqn{n \times n} matrix), \code{mu} (\eqn{\mu}, an
  \eqn{N}-element vector), \code{delta} (\eqn{\delta}, an
  \eqn{H}-element vector), \code{alpha}
  (\eqn{A}, an \eqn{N \times N} matrix),
  \code{var.add} (\eqn{\gamma}, an \eqn{H}-element vector),
  \code{var.mult} (\eqn{\nu}, an \eqn{H}-element vector), and
  \code{var.pow} (\eqn{\xi}, an \eqn{H}-element vector), where \eqn{N}
  is the number of series and \eqn{H} is the number of houses.
  Parameters can always be omitted from \code{init} and/or
  \code{restrict}, in which case reasonable guesses are taken.
  If \code{mu.by.process} is \code{TRUE}, then \code{mu} is a \eqn{K
    \times N} matrix, where \eqn{K} is the number of processes.
  If \code{delta.by.process} is \code{TRUE}, then \code{delta} should be a
  \eqn{K \times H} matrix.

  The actual estimation is carried out by the internal function
  \code{\link{monocar.cal}}. See \code{\link{monocar.cal}} for
  additional arguments which control the estimation process. The
  arguments \code{optimizers}, \code{maxiters},
  \code{tolerances}, \code{epsilon},
  \code{lower.bounds}, \code{upper.bounds}, and \code{var.centers}
  can be included as arguments to
  \code{monocar.estimate}, which passes them to \code{\link{monocar.cal}}.
}
\value{
  An object of class \code{monocar}.
}
\references{Tahk, A. M. (2015) \dQuote{A Continuous-Time,
    Latent-Variable Model of Time-Series Data,} \emph{Political
    Analysis} 23 (2): 278--298.}
\seealso{
  \code{\link{summary.monocar}}, \code{\link{monocar.estimate}},
  and \code{\link{create.ctdata}}.
}
\examples{
## Run a simple MONOCAR model
car.model <- monocar.estimate(data = PresUnemploy)
summary(car.model)
                                  
## Run a MONOCAR model with a Bush administration
##   dummy variable added
car.model.bush <-
    monocar.estimate(data = PresUnemploy, exo="Bush")
summary(car.model.bush)

## Run a MONOCAR model on Gallup polls only
car.model.gallup <-
    monocar.estimate(data = PresUnemploy,
                     subset = house == "Gallup" |
                         house == "unemployment")
summary(car.model.gallup)

## Run a MONOCAR model with correlated error processes
car.model.cor <-
    monocar.estimate(data = PresUnemploy,
                     restrict=list(sigma="unrestricted"))
summary(car.model.cor)

## Run a MONOCAR model with correlated error processes and transformed variances
car.model.cor.var <-
    monocar.estimate(data = PresUnemploy,
                     restrict=list(sigma="unrestricted",
                         var.add=c(TRUE,FALSE),
                         var.mult=c(FALSE,TRUE),
                         var.pow=c(FALSE,TRUE)))
summary(car.model.cor.var)
}
\keyword{monocar}
\keyword{models}
\keyword{ts}
