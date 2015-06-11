\name{simulateRFHC}
\alias{simulateRFHC}
\title{Simulate Finite Gibbs Hardcore Model with Smoothly Varying Range}
\description{
Simulate stationary and isotropic 2D Mat\'ern Gaussian random field using Markov approximation in a grid.

Then make exponential transformation, and get a field (say) R. 

Then simulate finite Gibbs hardcore model so that the repulsion radius at s is R(s). 
}
\usage{
simulateRFHC(n=100, iter=1e5, ncol=50, window=square(),	
             r.field.pars=c(mean=-2, tau=5, range=0.4, nu=2), 
             start.pattern=NULL, pre.r.field=NULL, giveup=1e6, 
             type="HC", dbg=0, seed, cyclic, expand)
}
\arguments{
  \item{n}{Target number of points to generate into the simulated pattern.}
  \item{iter}{Number of Metropolis-Hastings iterations.}
  \item{ncol}{Grid x-dimension. y-dimension is computed using \code{window} so that aspect ratio is 1.}
  \item{window}{Simulation window, in spatstat's \code{\link{owin}}-object format.}
  \item{pre.r.field}{Pre-computed R-field. In spatstat's \code{\link{im}}-format. If null we simulate it with given parameters using 
   the function \code{\link{simulateGMRF}}.}
  \item{giveup}{See \dQuote{Details}.}
  \item{type}{Type of repulsion. HC is hard-core, HS hard-spheres (not yet supported).}
  \item{seed}{Random seed. If not given a system time based random seed is used.}
  \item{start.pattern}{Starting pattern. If null we start with a single uniformly random point.}
  \item{dbg}{Verbose output during simulation.} 
  \item{r.field.pars, cyclic}{Parameters for the field simulation, see \code{\link{simulateGMRF}} and \dQuote{Details}.}
}
\details{
The function simulates a finite Gibbs hard-core model where the hard-core radius, say R, depends on the spatial location 
and varies smoothly in the window. We use an exponentially transformed Mat\'ern family Gaussian random field for the smooth variation. 

The (non-negative) R-field can be provided using parameter \code{pre.r.field}. It should be produced with 
the function \code{\link{simulateGMRF}}.  If a field is not provided, the members of vector \code{r.field.pars} 
are used to simulate the field: we first simulate GMRF using \code{\link{simulateGMRF}} to get (say) \eqn{Y}{Y}, and 
then transform this to R=exp(Y+mean), where \sQuote{mean} is given in \code{r.field.pars}.

Given the R-field, we simulate a finite Gibbs hard-core model where the 
hard-core radius at s is given by R(s) (with nearest grid point approximation). We start 
with an empty configuration x. Until x has n points in it, we

i1. Suggest a uniform random point x1 from the window.\cr
i2. If the union (x,x1) does not violate hardcore rule, accept the union as the new pattern.
 If violates (i.e. exists j s.t. ||x1-xj||<min(R(x1), R(xj)) ), discard x1.\cr
i3. If 1000 consecutive discards, remove a randomly chosen point from x.\cr
i4. Goto i1.

This loop is run until EITHER x has n points in it OR \code{giveup} amount of iterations has been run.

After the initial configuration is set, we move a single point, i.e. change location uniformly, and accept the move if hard-core rule is not violated. This is repeated \code{iter} times.

The first step (i1-i4) is a birth-and-death process simulation, and the second step is Metropolis-Hastings simulation, see reference. 
}
\value{
Object of class \code{rfhcSim} which is the same as object of class \code{\link{ppp}} with the additional elements

$parameters : the parameters provided\cr
$r.field : The R-field used in simulation.
}
\seealso{
\code{\link{plot.rfhcSim}}-method for plotting the field and the pattern simultaneously.
}
\references{
Rajala, T. and Penttinen, A. (2012): Bayesian analysis of a Gibbs hard-core point pattern model with varying
repulsion range, \emph{Comp. Stat. Data An.}

Moller, Waagepetersen (2004) \emph{Statistical Inference and Simulation for Spatial Point Processes} Chapman&Hall/CRC
}
\examples{
 y <- simulateRFHC()
 plot(y)
}
\keyword{RFHC}
     