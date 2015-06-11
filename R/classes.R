###################################################################
## different classes used by the package rfhc
##


rfhc_class_Sim<-function(x, ...){
	class(x)<-c("rfhcSim", is(x))
	x
}
rfhc_class_Fit<-function(x, ...){
	class(x)<-c("rfhcFit", is(x))
	x
}
rfhc_class_Field<-function(x, ...){
	class(x)<-c("rfhcField", is(x))
	x
}

rfhc_class_FitSummary<-function(x, ...){
	class(x)<-c("rfhcFitSummary", is(x))
	x
}
