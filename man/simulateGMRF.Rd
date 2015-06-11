\name{simulateGMRF}
\alias{simulateGMRF}
\alias{simulate_condGMRF}
\title{Simulate 2D Mat\'ern Field Using Markov Approximation.}
\description{
Simulate stationary and isotropic 2D Mat\'ern Gaussian random field using Markov approximation in a grid. Also conditional with given data.
}
\usage{
simulateGMRF(ncol=50, nrow=50, tau=3, range=0.3, nu=2, expand,
             xlim=c(0,1), ylim=c(0,1), cyclic=FALSE, dbg=0)
simulate_condGMRF(ncol=50, nrow=50, tau=5, range=0.2, nu=2, 
                  expand, xlim=c(0,1), ylim=c(0,1), 
                  data=list(x=.5, y=.5, v=0), 
                  cyclic=FALSE, dbg=0, use.RandomFields=FALSE, ...)
}
\arguments{
  \item{ncol, nrow}{Grid dimensions. ncol in x direction.}
  \item{xlim, ylim}{Spatial dimensions, rectangular window.}
  \item{expand}{Expansion of the grid for border bias elimination. If not given, computed using \code{nu}.}
  \item{tau}{Precision i.e. inverse variance.}
  \item{range}{Range of covariance, controls smoothness.}
  \item{nu}{The differentiability parameter of covariance. Accepts: 1,2,3.}
  \item{cyclic}{Use cyclic (toroidal) grid.}
  \item{dbg}{Verbose output during simulation.}
  \item{use.RandomFields}{Option to use randomFields package for conditional simulation.}
  \item{...}{If using randomFields, these are passed on to the CondSimu-function.}
  \item{data}{List with components x,y and v giving the x- and y-coordinates of values v for conditional simulation.}  
}
\details{
The stationary and isotropic Mat\'ern covariance function at distance t is

\deqn{C(t) \propto tau(kappa t)^nu K_{nu}(kappa t)}{C(t) \propto tau*(kappa*t)^nu K_{nu}(kappa*t)}

where  \eqn{\tau}{tau} is the precision (inverse variance), \eqn{\nu}{nu} is differentiability of the covariance, 
and \eqn{\kappa}{kappa} is the scale. We use a transformed kappa called range \eqn{r=\frac{\sqrt{8\nu}}{\kappa}}{r=sqrt(8*nu)/kappa}
so that \eqn{C(r)\approx 0.1}{C(r)~=0.1}. 

The simulation is done using Markov approximation, see Held \& Rue 2005 and Lindgren et al 2011. Hence 
the term GMRF: Gaussian Markov Random Field.
}

\value{
Object of class \code{rfhcField} which is the same as object of class \code{\link{im}} with the simulation parameters 
stored in the element $parameters.
}
\references{
Held, Rue (2005) \emph{Gaussian Markov Random Fields: Theory and Applications}.
       Chapman&Hall/CRC.

Lindgren, Rue, Lindst\"om (2011) \emph{An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach}
	J. R. Statistic. Soc. B. 73.	
}
\examples{
 y<-simulateGMRF()
 plot(y)
}
\keyword{GMRF}
     