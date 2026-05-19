\name{summary.monocar}
\alias{print.monocar}
\alias{summary.monocar}
\alias{print.summary.monocar}
\title{Print and summarize MONOCAR models}
\description{Print and summary methods for MONOCAR
  models from \code{monocar.estimate}.}
\usage{
\method{print}{monocar}(x, reverse.offdiag = TRUE,
    digits = max(3L, getOption("digits") - 3L), \dots)

\method{summary}{monocar}(object, reverse.offdiag = TRUE, \dots)

\method{print}{summary.monocar}(x,
     digits = max(3L, getOption("digits") - 3L),
     signif.stars = getOption("show.signif.stars"), \dots)
}
\arguments{
  \item{x}{For \code{print.monocar}, an object of class \code{monocar}.
    For \code{print.summary.monocar}, an object of class
    \code{summary.monocar}.}
  \item{object}{An object of class \code{monocar}.}
  \item{reverse.offdiag}{A logical value indicating whether to reverse
     the sign of off-diagonal elements of Theta. In a MONOCAR model, a
     positive value of an off-diagonal element of Theta indicates that
     an increase in one latent process tends to \emph{decrease} the
     level of another latent process. Reversing the sign of these
     parameters means that a positive value indicates one process tends
     to increase another. This can simplify the interpretation of the
     parameters. The default is to reverse the sign of off-diagonal
     elements of Theta. For \code{print.summary.monocar}, whether these
     elements are reversed is set by \code{summary.monocar} and stored
     in the \code{summary.monocar} object.}
  \item{digits}{The minimal number of significant digits, see
  \code{print.default}.}
  \item{signif.stars}{A logical value indicating whether stars should be
    printed based on p-values, see \code{printCoefmat}.}
  \item{\dots}{Additional arguments to pass to \code{print.default} or
    \code{printCoefmat}.}
}
\value{
  \code{summary.monocar} returns an object of class
     \code{summary.monocar} containing a summary information about the
     model.
}
\seealso{
  \code{\link{monocar.estimate}} and \code{\link{monocar.hist}}.
}
\examples{
## Run a simple MONOCAR model
car.model <- monocar.estimate(data = PresUnemploy)
## Print model estimates
car.model
## Print model estimates with diagonal not reversed
print(car.model, reverse.offdiag=FALSE)
## Print summary information about the model
summary(car.model)
}
\keyword{monocar}
\keyword{ts}
