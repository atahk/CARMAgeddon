\name{monocar.cal}
\alias{monocar.cal}
\title{Internal routine to provide estimate of MONOCAR model}
\description{
  Internal routine to provide estimate of a Multivariate overlapping
  noisy observation continuous-time autoregressive (MONOCAR) model
}
\usage{
monocar.cal(init.param, data, exodata, pc.list, verbose=1,
  optimizers = c("MMA", "BFGS"), maxiters = NULL, tolerances = NULL,
  epsilon = 0, lower.bounds = NULL, upper.bounds = NULL,
  var.centers = NULL, priors = NULL, newton.raphson = 0,
  time.indices = NULL, compute.hessian = TRUE,
  compute.partial.scores = TRUE, history.times = NULL,
  use.parallel = NULL)
}
\arguments{
  \item{init.param}{A vector giving the initial parameter guesses.}
  \item{data}{A data frame giving the data (excluding exogenous
    variables).}
  \item{exodata}{A data frame giving the data for exogenous variables
    (if any are present).}
  \item{pc.list}{A list of matrices, one for each parameter, giving the
    parameter restrictions.}
  \item{verbose}{An integer indicating the verbosity of the
    output. Larger numbers mean greater verbosity. Zero means no output
    while running.}
  \item{optimizers}{A character vector giving the names of optimization
    algorithms, in the order in which
    they should be run. See details for more information.}
  \item{maxiters}{The maximum number of iterations for each optimizer
    (in the same order as \code{optimizers}). A negative value means no
    maximum.}
  \item{tolerances}{The tolerance used for each optimization algorithm
    (in the same order as \code{optimizers}).}
  \item{epsilon}{The value of \eqn{\epsilon} to be used in calculating
    numerical derivatives and Hessians. If equal to zero, analytic
    derivatives and Hessians are used instead. The default value is
    zero, meaning analytic derivatives and Hessians. This is usually
    preferable to a numerical approach.}
  \item{lower.bounds}{A vector giving the lower bound for each element
    of \code{init.param}}
  \item{upper.bounds}{A vector giving the upper bound for each element
    of \code{init.param}}
  \item{var.centers}{For each house, the value around which to recenter
    the variances. See \code{\link{monocar.estimate}} for more
    details. This is the \eqn{\bar{v}} that appears in the
    observation equations.}
  \item{priors}{An optional matrix specifying prior distributions on
    beta for Bayesian maximum a posteriori estimation.}
  \item{newton.raphson}{Number of Newton-Raphson interations using the
    analytic Hessian as an alterative to or in additon to optimizers
    from nloptr. Performed after other optimizers. The default is zero.}
  \item{time.indices}{An optional vector specifying time indices.}
  \item{compute.hessian}{Whether to compute the Hessian. This is needed
    to compute the variance-covariance matrix and standard errors for
    the parameters. The default is to compute the Hessian.}
  \item{compute.partial.scores}{Whether to compute the partial
    derivative of the log-likelihood contribution of each observation
    with respect to the parameters. The default is to compute these
    quantities.}
  \item{history.times}{An optional vector of times at which to compute the
    estimated history of the latent processes. The default of \code{NULL}
    means the history is not computed.}
  \item{use.parallel}{Whether to compute the gradient using parallel
    processing. This is faster than the non-parallel machines on
    multi-core machines and the same speed on single-core
    machines. However, computational time is not linear in the number of
    cores. If estimating several models simultaneously, it may be faster to
    dedicate a single core to each, set \code{use.parallel} to
    \code{FALSE}, and estimate the models in parallel. The default is
    \code{NULL}, which is interpreted as \code{TRUE} if and only if
    the number of threads is greater than one and \code{FALSE}
    otherwise.}
}
\details{
  This is an internal function used to estimate MONOCAR models.
  Parameters can be passed to it from \code{\link{monocar.estimate}} or
  \code{\link{monocar.hist}}.
  It is normally only called internally by
  \code{\link{monocar.estimate}} or
  \code{\link{monocar.hist}} and not by the user.

  Possible options for \code{optimizers} are:
  \code{MMA} (Method of Moving Asymptotes),
  \code{BFGS} (low-storage Broyden-Fletcher-Goldfarb-Shanno),
  \code{NM} (Nelder-Mead),
  \code{DIRECTL} (DIviding RECTangles, locally biased),
  \code{CRS} (Controlled Random Search with local mutation),
  \code{MLSL} (Multi-Level Single-Linkage),
  \code{STOGO} (StoGO),
  \code{STOGORAND} (randomized variant of StoGO),
  \code{ISRES} (Improved Stochastic Ranking Evolution Strategy),
  \code{ESCH} (evolutionary algorithm),
  \code{COBYLA} (Constrained Optimization BY Linear Approximations),
  \code{BOBYQA} (BOBYQA),
  \code{NEWUOA} (NEWUOA),
  \code{PRAXIS} (PRincipal AXIS),
  \code{SBPLX} (Sbplx, based on Subplex),
  \code{SLSQP} (Sequential Least-Squares Quadratic Programming),
  \code{TNEWTON} (Preconditioned truncated Newton),
  \code{VAR2} (Shifted limited-memory variable-metric, rank 2),
  and
  \code{VAR1} (Shifted limited-memory variable-metric, rank 1).
  The default of \code{MMA} is usually provides a good mix of reliability and speed.
  \code{BFGS} is faster, but does not always behave as well when the
  initial guess is far from the maximum. \code{NM} is slower than either
  of these, but can be more reliable when the initial guess is far from
  the maximum. Beginning with more reliable algorithms followed by
  faster algorithms can be useful in difficult cases.  
  See the
  \href{http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms}{NLopt
  documentation} for details on these algorithms.
}
\value{
  A list containing the estimated parameters, Hessian, and other information.
}
\references{Tahk, A. M. (2015) \dQuote{A Continuous-Time,
    Latent-Variable Model of Time-Series Data,} \emph{Political
    Analysis} 23 (2): 278--298.}
\seealso{
  \code{\link{monocar.estimate}} and \code{\link{monocar.hist}}
}
\keyword{internal}
\keyword{monocar}
\keyword{models}
\keyword{ts}
