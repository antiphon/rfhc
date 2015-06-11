\name{methods_rfhcFitSummary}
\alias{plot.rfhcFitSummary}
\alias{print.rfhcFitSummary}
\alias{summary.rfhcFit}
\alias{rfhcFitSummary}
\title{Plot and Print rfhcFitSummary-objects.}
\description{
Plot and print the result object of function \code{\link{summary.rfhcFit}}.
}
\usage{
\method{summary}{rfhcFit}(x, burnin=1, thin=1, Rburnin=1, Rthin=1, paradj=1, Radj=1, ...)

\method{print}{rfhcFitSummary}(x, ...)
\method{plot}{rfhcFitSummary}(x, ..., maxcols=4)
}
\arguments{
  \item{x}{in \sQuote{summary}: Result of \code{\link{fitRFHC}}, see \code{\link{rfhcFit}}. Otherwise the result of the \sQuote{summary}-function.}
  
  \item{maxcols}{Adjust plot layout. Maximum number of columns.}  
  \item{burnin, thin, Rburnin, Rthin, paradj, Radj}{See \code{\link{rfhcFit}}.}
  \item{...}{In plot : Additional parameteres passed on to the function \code{\link{plot.im}} that plots the field. E.g. colors, ribbon options. Others: not used.}
}
\details{
	Write R=Robs union Rs, where Robs is the dataholding R's and Rs is predicted outside data.

	Print: Prints some useful info.

	Summary object contains the following general elements:
	
	$parameters: the hyperparameters and beta given for the fitting\cr
	$options: options used in the fitting\cr
	$sampling: sampling scheme (burn-in, thinning)\cr
	$estimated: names of estimated hyper parameters.
	
	For each estimated hyperparameter, say tau, the following sub-elements are present:
	
	$summary : basic mean/sd/quantile table\cr
	$density : estimated smoothed posterior density using the function \code{\link{density}} 
	
	$Robs :  contains the $summary table and $density which is a list of estimated smoothed posterior densities.
	
	$Rpred : $mean.im and $sd.im are \code{\link{im}}-objects for predictive mean and sd of the latent R. $var has $summary and $density, 
	for the pointwise variances of log-R values.  
	
	
	
	Plot: plots the posteriors of estimated hyper-parameters and the posterior mean and pointwise-variance of R-field. 
	
}
\seealso{
\code{\link{fitRFHC}}, \code{\link{rfhcFit}}.
}
\keyword{rfhc}
     