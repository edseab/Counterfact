
counterfact <- function (object, x, range = "standard", n=1000, default = mean, CI=FALSE, PI=FALSE, other=NULL, CIsims = 1000, PIsims= 1000, conf=95, odds.to.prob=F, data=summary(object)$call$data, progress="none",...){


progress <- function (i,length, percent=1) {
	if(i==1){print(0)}
	if(trunc(i*100/length) %% percent ==0 & trunc(i*100/length)!=trunc((i-1)*100/length)){ print(trunc(i*100/length))}
}

# Classic errors
if(is.null(data)) stop("Model data must be inputed as data frame")
if(!is.character(x)) stop("x must be inputed as character object")

d <- get(as.character(data))
# if (any(apply(d,2,function(x) !is.numeric(x)))){
	# stop("Non numeric variable - choose counterfactual value")}

# extract range	
if (length(range)==1){if(range=="standard") range <- range(d[,x],na.rm=T)}
f <- as.function (default)

ifelse(class(object)[1]=="lmerMod" | class(object)[1]=="glmmPQL",
				betas <- fixef(object),
				betas <- coef(object))
vars <- attr(terms(object),"term.labels")
coefs<-vars
if(length(grep(":",vars))>0) coefs <- vars[-grep(":",vars)]
if(length(grep("I\\(",coefs))>0) coefs <- coefs[-grep("I\\(",vars)]
d <- d [,coefs]

factors <- names(which(unlist(lapply(d,function(.)is.character(.)|is.factor(.)))))
if (length(factors)>0){
	if(any(!(factors %in% names(other)))) 
					stop (paste("Non numeric variable: (",paste(factors,collapse=";"),"). Choose counterfactual value"), sep="")
	}

nd <- unlist(lapply(d[,which(!(colnames(d) %in% factors))],function(x)f(x,na.rm=T)))
if(length(factors)>0) {for (i in 1:length(factors)){
	nd <- append(nd, NA, which(colnames(d)==factors[i])-1)}
	names(nd)[which(colnames(d) %in% factors)]<-factors}


if (!is.null(other)) {
	if (!is.list(other)) stop ("other variables must be inputed as list")
	if (any(!(names(other) %in% coefs))) {
		wrong <- names(other)[which(!(names(other) %in% coefs))]
		if (length(wrong)>1) stop (paste("Variables", wrong, "not found in model"))
		if (length(wrong)==1) stop (paste("Variable", wrong, "not found in model"))
		}
	ph <- other[which(names(other) %in% coefs)]
	nd [match(names(ph), coefs)] <- as.vector(unlist(ph))
	}
nd <- data.frame(matrix(nd,n,length(nd),byrow=T))
nd[,which(!(colnames(d) %in% factors))] <- unlist(lapply(as.list(nd)[which(!(colnames(d) %in% factors))],function(.)as.numeric(as.character(.))))
nd[,which(colnames(d) %in% factors)] <- unlist(lapply(data.frame(nd[which(colnames(d) %in% factors)]),as.character))
colnames (nd) <- coefs

xvalues <-seq(range[1],range[2],length.out=n)
nd [,grep(x,coefs)] <- xvalues

if(length(factors)>0){
	extras <- prod(unlist(lapply(as.list(d)[factors],function(.)length(unique(.)))))
	nd <- rbind (nd, nd[(n-extras+1):n,])
	lvls <- lapply(as.list(d)[factors],function(.)levels(as.factor(.)))
	ifelse (length(lvls)>1,
		for(i in 1:length(lvls)) nd[(n+1):(n + extras),names(lvls)[i]] <-as.vector(unlist(lvls[i])),
		nd[(n+1):(n + extras),names(lvls)] <- as.vector(unlist(lvls)))
}
	
nd <-   as.list(nd)


mm <- model.matrix(as.formula(paste("~",paste(vars, collapse="+"))), data=nd)
mm <- mm[1:n,]
y_hat<-mm%*%betas

pr <- data.frame(cbind(xvalues,y_hat))
colnames (pr) [2] <- "predicted.mean"

if (CI) {
if(class(object)[1]=="lmerMod"){ predFun<-function(.) mm%*%fixef(.)
bb<-bootMer(object,FUN=predFun,nsim=CIsims, .progress=progress)
bb_se<-apply(bb$t,2,function(x) x[order(x)])
lo <- which.min(abs(1:CIsims - (100-conf)*CIsims/200))
hi <- CIsims - lo
pr$LowerCI <- bb_se[lo,]
pr$UpperCI <- bb_se[hi,]
}

if(class(object)[1]=="lm"){
v <-vcov(mod)
var.pred <- rowSums((mm %*% v) * mm) # more efficient way of calculating diag(mm %*% v %*% t(mm))
# this is because var(y_hat) = var(X*B_hat) = X*var(B_hat)*t(X)
se.pred<-sqrt(var.pred)
tval <- qt((100-conf)/200, df=(df.residual(mod)))
pr$LowerCI <- y_hat+se.pred*tval
pr$UpperCI <- y_hat-se.pred*tval
}

if(class(object)[1]=="glmmPQL"){
	vc <- vcov(object)
	sim <- matrix(NA, CIsims, n)
	for (i in 1: CIsims){
		betas.sim <- rmvnorm(n=1,mean=betas,sigma=vc)
		sim[i,] <- mm%*%as.vector(betas.sim)
	}
	pr$LowerCI <- apply(sim,2,quantile,(1-conf/100)/2)
    pr$UpperCI <- apply(sim,2,quantile,1-((1-conf/100)/2))
}
}
return (pr)


}
