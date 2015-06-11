# Methods for class rfhcSim (pp simulation)
# 
# Author: antiphon
###############################################################################


plot.rfhcSim<-function(x, size=0, field=TRUE, ...) {
	if(field)
		plot(x$r.field, ...)
	else plot.ppp(x, ...)
	points(x$x, x$y, cex=0.3, pch=19, col="goldenrod")
	if(size>0) symbols(x$x, x$y, circles=x$marks*size, inches=FALSE, add=TRUE)
}

print.rfhcSim<-function(x, ...) {
	cat("RFHC-simulation\n")
	print.ppp(x, ...)
	cat("r-field parameters:", paste(names(x$parameters$r.field.pars),x$parameters$r.field.pars, sep=":"),"\n")
	cat("iterations:", x$parameters$iterations,"\n")
}