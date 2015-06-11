\name{methods_rfhcSim}
\alias{plot.rfhcSim}
\alias{print.rfhcSim}
\alias{rfhcSim}
\title{Plot and Print rfhcSim-objects.}
\description{
Plot and print the result object of function \code{\link{simulateRFHC}}.
}
\usage{
\method{plot}{rfhcSim}(x, size=0, field=TRUE , ...)
\method{print}{rfhcSim}(x, ...)
}
\arguments{
  \item{x}{Simulated pattern.}
  \item{size}{If >0, draw the radius for each point multiplied by \sQuote{size}. Best use: size=1.}
  \item{...}{in plot: Additional parameteres passed on to the function \code{\link{plot.im}} that plots the field. E.g. colors, ribbon options.}
}
\details{
	Plot: First plot the field using \code{\link{plot.im}}, then overlay the points using \code{\link{points}}. If size=1 use \code{\link{symbols}} to add circles.
	
	Print: Print short description.
}
\seealso{
\code{\link{simulateRFHC}}.
}
\examples{
 y <- simulateRFHC()
 \dontrun{plot(y, 1, col=gray.colors(30))}
 print(y)
}
\keyword{rfhc}
     