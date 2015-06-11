\name{fitRFHC}
\alias{fitRFHC}
\title{Fit Finite Gibbs Hardcore Model with Smoothly Varying Range}
\description{
Finite Gibbs hardcore model has repulsion radius at s given by some latent field R(s). Given a point pattern x, we estimate R using Bayesian MCMC. 
}

\usage{
fitRFHC(x, iter=100, ncol=50, 
		nu=2,
		range=0.3,
		mupars=c(free=1, init=NA, m=0, s2=100),
		taupars=c(free=1, init=5, F=1.05, alpha=0.25, beta=0.005),
		Rpars=list(choose=0.15, from=NULL), 
		beta=9,  
		expand, cyclic,
		seed, dbg=1, RtoNND=0.999, print.prefix="", type="HC",
		predict.R=10, 
		predict.x=0)
}
\arguments{
  \item{x}{Data point pattern, best given in the spatstat's \code{\link{ppp}}-object format.}
  \item{iter}{Number of iterations in the MCMC.}
  \item{ncol}{Latent field grid x-dimension. y-dimension is computed using \code{window} object of \code{x} so that aspect ratio is 1.}
  \item{nu}{Differentiability (smoothness) parameter of Mat\'ern covariance. Supports 1,2,3.}
  \item{range}{The range parameter.}
  \item{mupars}{Vector of parameters for the mean parameter, see \dQuote{Details}.}
  \item{taupars}{Vector of parameters for the precision parameter, see \dQuote{Details}.}
  \item{Rpars}{Vector of parameters modifying the fitting of the latent field, see \dQuote{Details}.}
  \item{beta}{The expected packing density parameter, see \dQuote{Details}.}
  \item{predict.R}{Frequency for predicting R-values outside datapoints using conditional simulation.}
  \item{predict.x}{Frequency for predicting the point process, see \dQuote{Details}.}
  \item{RtoNND}{Initial R(s), s in x, values are nearest neighbour distances multiplied with \sQuote{RtoNND}.}
  \item{type}{Type of repulsion. HC is hard-core, HS is hard-spheres (not yet supported).}
  \item{seed}{Random seed. If not given a system time based random seed is used.}
  \item{dbg}{Verbose output during MCMC. From 0 upwards.} 
  \item{print.prefix}{Prefix for run-time messages when dbg>0. (useful for large estimation runs with several nu, beta, etc. parameters).}
  \item{expand, cyclic}{Parameters for the field simulation, see \code{\link{simulateGMRF}}.}
}
\details{
The function fits a finite Gibbs hard-core model where the hard-core radius, say R, depends on the spatial location and varies smoothly in the window. We use 
an exponentially transformed Mat\'ern field as the prior for the smooth variation. 

Give data in \code{\link{ppp}}-format (or as a result of \code{\link{simulateRFHC}}). The size of data is denoted by \sQuote{n} in what follows. 

We estimate the latent, smoothly varying radius field R using a log-Gaussian prior. The Gibbs likelihood is replaced with a \sQuote{modified pseudo-likelihood}. 
Instead of using the Papangelou-intensity, which involves integration over the unknown R, we model the discrepancy between the nearest neighbour distances of data
and the R-values at data locations. If \sQuote{Ri} denotes the R(x_i) value and \sQuote{di} the nearest neighbour distance of point \sQuote{x_i},
the connection between the R and d values is modelled assuming independence of the terms
\deqn{\frac{d_i-R_i}{d_i}\sim Beta(1, beta) .}{(di-Ri)/di ~ Beta(1, beta) .} 

The \sQuote{beta} parameter controls the amount of discrepancy we allow between \sQuote{Ri} and \sQuote{di}.  If set to 1 the model only says that 'Ri is between 0 and di', and if set very high
we assume that 'Ri is di'.   In practice it should be between 5-25, amounting to approximate packing densities, i.e. the amount of points the field can support, of 50-100\%. The number depends also on the smoothness of the field.
See the referenced article for more details.

The hyper priors are given in vector form as follows:

* mupars=c(free=1, init=NA, m=0, s2=100): First element (\sQuote{free}) flags whether the mean parameter is to be estimated. The initial value (\sQuote{init}) can be missing 
and if so, is set to the mean of the initial log R-values (see below). The prior for mean is Gaussian with parameters \sQuote{m} and \sQuote{s2}. We use a Gibbs step so no need to
define proposal distribution. 

* taupars=c(free=0, init=5, F=1.05, alpha=0.25, beta=0.005): First two like in mupars. For proposal we use \eqn{\tau'=f \tau}{tau'=f*tau} where \sQuote{f} is Lambert distibuted with parameter \sQuote{F}.
The prior is Gamma(alpha, beta). 

Note 1: Fixed values can be given also in the form \sQuote{mupars=c(free=0, <init>)} or just \sQuote{mupars=<init>} where \sQuote{<init>} is your fixed value. Same for tau.

Note 2: The range is best set to some fixed value as data has very little info on it. If you must, you can make it free using the syntax of \sQuote{taupars}. Be aware that this increases
computations enormously as we need to factorise the n^2 precision matrix each step.

We primarily estimate the R values at data points, say \sQuote{Robs=(R1, R2,...,Rn)}. Then when a prediction outside datapoints, say \sQuote{Rs} with union(Rs, Robs)=R, 
is requested every min(predict.x, predict.R) step, we do a conditional simulation to \sQuote{Rs} given \sQuote{Robs}.

We also sample from the marginal posterior of the point process x every \sQuote{predict.x} step. When this is requested, the algorithm first forces prediction of \sQuote{Rs}, 
and then uses \code{\link{simulateRFHC}} with the current \sQuote{R=(Robs, Rs)} as parameter \sQuote{pre.r.field}. Point count n and window are like in data.

The initial \sQuote{Robs} values are set to nearest neighbour distances multiplied with \sQuote{RtoNND}. 

Each iteration we attempt to update \sQuote{Rpars$choose} number of randomly chosen \sQuote{Robs} points. If \sQuote{Rpars$from} is NULL it is taken as the vector (1,...,n). Otherwise you can give a 
subset to keep some locations fixed to their initial values.

}
\value{
Object of class \code{\link{rfhcFit}}. It contains the following elements:

$mu, $tau : The histories. If fixed, single value.\cr
$range : Fixed value. If free, the history.\cr
$Robs : (iter x n) matrix. Each row is the state of Robs, so each column i is the history of R(x_i).\cr
$Rpred : Matrix. Each row is the prediction Rs union Robs. Number of rows depends on predict.R and predict.x.\cr
$time : Time used.\cr
$pphist : List with x predictions.\cr
$options: Options set for the algorithm excluding hyper parameters, R update and beta.\cr
$dataholders : Indices of R-grid (in columnwise order) of those R-values that hold a datapoint, i.e. R[dataholders]=Robs.\cr
$acceptance : Number of accepted steps, separate for Robs and hyper parameters.\cr
$parameters : The hyper parameter options, R update options and beta.\cr 
$data : The given data object, \sQuote{x}.\cr
$Rlast : Last state of R as an \code{\link{im}}-object. Note that depending on the \sQuote{predict.R/predict.x} parameter the Rs values might not be up-to-date with Robs. 
}

\seealso{
\code{\link{plot.rfhcFit}}-method, \code{\link{summary.rfhcFit}}-method, \code{\link{print.rfhcFit}}-method, \code{\link{plotRtraj.rfhcFit}}-function, \code{\link{plotRpost.rfhcFit}}-method.
}
\references{
Rajala, T. and Penttinen, A. (2012): Bayesian analysis of a Gibbs hard-core point pattern model with varying
repulsion range, \emph{Comp. Stat. Data An.}
}

\examples{
 x <- simulateRFHC()
 y <- fitRFHC(x)
 \dontrun{plot(y)}
}
\keyword{rfhc}
     