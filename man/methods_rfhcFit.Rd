\name{methods_rfhcFit}
\alias{plot.rfhcFit}
\alias{print.rfhcFit}
\alias{plotRtraj.rfhcFit}
\alias{plotRpost.rfhcFit}
\alias{rfhcFit}
\title{Plot and Print rfhcFit-objects.}
\description{
Plot and print the result object of function \code{\link{fitRFHC}}.
}
\usage{
\method{plot}{rfhcFit}(x, burnin=1, thin=1, Rburnin=1, Rthin=1, data=TRUE, sigma=FALSE, ...)
\method{summary}{rfhcFit}(x, burnin=1, thin=1, Rburnin=1, Rthin=1, paradj=1, Radj=1, ...)

\method{print}{rfhcFit}(x, ...)

plotRtraj.rfhcFit(x, these, burnin=1, thin=1, ...)
plotRpost.rfhcFit(x, burnin=1, thin=1, these=1:6, maxcol=6, breaks=50, error=TRUE, ...)
}

\arguments{
  \item{x}{Result of \code{\link{fitRFHC}}.}
  \item{burnin, thin}{The posterior esimates are based on a MCMC sample taken with vector \code{seq(burnin, iter, by=thin)}.}
  \item{Rburnin, Rthin}{Burn-in and thinning for the predicted R-values, i.e. we sample the prediction chain with \code{seq(Rburnin, Riter, by=Rthin)} where \sQuote{Riter} depends on \sQuote{predict.R} (see \code{\link{fitRFHC}}).}
  \item{data}{Should the data be plotted.}
  \item{sigma}{If tau estimated, plot the history as tau or sigma^2=1/tau.}
  \item{paradj, Radj}{Bandwith adjusting parameter given to \code{\link{density}} when computing posterior densities, see \dQuote{Details}. \sQuote{paradj} for hyperparameters, \sQuote{Radj} for R(x_i)-values.}
  \item{these}{Which of the dataholding R values should be plotted. Indexing matches that of data points i.e. i=1,...,n.}
  \item{maxcol, breaks}{Layout width and resolution of histogram.}
  \item{error}{If true, plot the \eqn{(d_i-R_i)/d_i}{(di-Ri)/di} instead of \eqn{R_i}{Ri}.} 
  \item{...}{Additional parameteres passed on to the function \code{\link{plot.im}} that plots the field. E.g. colors, ribbon options.}
}

\details{
	Write R=Robs union Rs, where Robs is the dataholding R's and Rs is predicted outside data.
	
	Print: Print short description. Note especially the length of R chain, this is important for choosing \sQuote{Rburnin} and \sQuote{Rthin}.
		
	Plot: Plots the data, the MCMC histories of estimated hyperparameters and the posterior mean field + the median of posterior Robs values. If the library \sQuote{mmap} is installed (from the author),
	the median marks are compared to true values and a kernel smoothed mean field of the squared difference is plotted.
	
	Plot R trajectories: Plots the history of Robs values.
	
	Plot R posterior: Plots the posterior histogram of Robs values. Either the value or the relative error to nearest neighbour distance. 
	Red line is the beta-link (see \code{\link{fitRFHC}}).
	
	Summary: Summarises the fit. Value is an object of class \code{\link{rfhcFitSummary}}, see the class-page for more details.
	
}
\seealso{
\code{\link{fitRFHC}}, \code{\link{rfhcFitSummary}}.
}
\examples{
 x <- simulateRFHC()
 f <- fitRFHC(x)
 s <- summary(f)
 \dontrun{plot(s)}
 print(s) 
}
\keyword{rfhc}
     