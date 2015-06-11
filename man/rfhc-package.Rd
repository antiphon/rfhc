\name{rfhc-package}
\alias{rfhc-package}
\alias{rfhc}
\docType{package}
\title{
Finite Gibbs Hardcore Model with Smoothly Varying Range
}
\description{
Simulate and fit finite Gibbs hardcore model with local repulsion radii from a latent random field.
}
\details{
Most important functions are

\code{\link{simulateRFHC}}

for simulating the model,

\code{\link{fitRFHC}}

for fitting the model.

The dataset used in our example is included, see \code{\link{riparia}}.
}
\author{
Tuomas Rajala

Maintainer: Tuomas Rajala <tuomas.rajala@iki.fi>
}
\references{
Rajala, T. and Penttinen, A. (2012): Bayesian analysis of a Gibbs hard-core point pattern model with varying
repulsion range, \emph{Comp. Stat. Data An.}
}
\keyword{ package }
\seealso{
\code{\link[spatstat:spatstat-package]{spatstat}}

\code{\link[randomFields:randomFields-package]{randomFields}}
}
\examples{
x <- simulateRFHC(dbg=1)
f <- fitRFHC(x, dbg=1)
\dontrun{plot(f)}
s <- summary(f, burnin=40, thin=2)
print(s)
}
