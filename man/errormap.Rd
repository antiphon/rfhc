\name{errormapfun}
\alias{error.mapfun}
\alias{error.sumfun}
\title{rfhc: Compute IMSE from predictive simulation}
\description{
If the \code{\link{fitRFHC}} had predict.x>0, we can use the predictions to 
compute a model diagnostic.
}

\usage{
error.mapfun(x, burnin=1, thin=1, mapfun=distmap, ...) 
}

\arguments{
  \item{x}{Result of \code{\link{fitRFHC}}.}
  \item{burnin, thin}{How the predictive simulations of the point pattern are sampled. Note that the number of samples is dependent on \code{predict.x}.}
  \item{mapfun}{The function that maps a point pattern to a field. Should result in a \code{\link{im}}-object (we use the im-objects $v-matrix). 
    See e.g. \code{\link{distmap}} and \code{\link{density.ppp}} in package spatstat.}
  \item{...}{Additional parameters passed to \sQuote{mapfun}.}
}

\details{
	Let x0 be the data point pattern and x1,...,xk the predicted patterns.
	
	We compute mapfun(x0) to get v0 field, and compare this to v1,...,vk from the predictions.
	
	Summary indices computed from the fields are:
	\deqn{IMSE=a\sum_i MSE_i}{IMSE=a*sum(MSE_i)} \deqn{MSE_i=(1/k)\sum_j (v0_i - vj_i)^2}{MSE_i=(1/k)sum((v0_j - vj_i)^2)}
	where a is the area of a grid cell. 
	
	In IMSE_std the MSE_i is divided by the sd of vj_i.
	
	The result list contains also the maps used in computation.
	
	We discuss only the $IMSE value in the paper. The rest of the results can be ignored at this point.
	
	The idea is that this function is run for different fits with different beta and range parameters, and the IMSE values compared. 
}
\seealso{
\code{\link{fitRFHC}}, \code{\link{rfhcFit}}.
}

\references{
Rajala, T. and Penttinen, A. (2012): Bayesian analysis of a Gibbs hard-core point pattern model with varying
repulsion range, \emph{Comp. Stat. Data An.}
}


\examples{
 x<-simulateRFHC()
 f<-fitRFHC(x, predict.x=10)
 err <- error.mapfun(f, mapfun=density)
 err$IMSE
 }
 
\keyword{rfhc}