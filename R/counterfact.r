
counterfact <- function (object, x, range = "standard", n=1000, values = NULL, default = mean, CI=FALSE, PI=FALSE,
other=NULL, CIsims = 1000, PIsims= 1000, conf=95, unlink=F, data=summary(object)$call$data, progress="none", use.predict=F,mixture=FALSE,...){

if(mixture & class(object)[1]!="glmmTMB"){
  stop("Mixture models (zero-inflation, hurdle...) currently only compatible with glmmTMB models")
}
ifelse(class(object)[1]=="brmsfit", d <- object$data, d <- eval(parse(text=as.character(c(data)))))
if (!is.character(x)) x <- deparse(substitute(x))

# Classic errors
if (is.null(data) & class(object)[1]!="brmsfit") stop("Model data must be inputed as data frame")
if (!is.null(other) & !is.list(other)) stop("other variables must be inputed as list")

# extract range	
if (length(range)==1){if(range=="standard") range <- range(as.numeric(d[,x]),na.rm=T)}
f <- as.function (default)

if(class(object)[1] %in% c("brmsfit")){
 betas <- as.data.frame(fixef(object))$Estimate
 coefs <- names(object$data)[!names(object$data) %in% c(object$formula$resp, object$ranef$group, "Intercept")]
 vars <- attr(terms(as.formula(object$formula)), "term.labels")
 vars <- vars[!grepl(" | ", vars)]
}else {
if(class(object)[1] %in% c("lmerMod","glmmPQL","glmerMod","lmerModLmerTest")) betas <- fixef(object)
if(class(object)[1] %in% c("lm","glm")) betas <- coef(object)

vars <- attr(terms(object),"term.labels")
coefs<-vars
if(length(grep(":",vars))>0) coefs <- vars[-grep(":",vars)]
if(length(grep("I\\(",coefs))>0) coefs <- coefs[-grep("I\\(",vars)]
}
d <- d [,coefs]

factors <- names(which(lapply(d,class) == "character" | lapply(d,class) =="factor"))
if (length(factors)>0){
    unspecced  <- !(factors %in% names(other))
	if(any(unspecced)) 
					stop (paste("Non numeric variable: (",paste(factors[unspecced],collapse=";"),"). Choose counterfactual value"), sep="")
	}

if(all(colnames(d) %in% factors)) stop ("counterfact bugs with no continuous predictors. For now.")

nd <- unlist(lapply(as.data.frame(d[,which(!(colnames(d) %in% factors))]),function(x)f(x,na.rm=T)))
names(nd) <- colnames(d)[which(!(colnames(d) %in% factors))]
if (length(factors)>0) {for (i in 1:length(factors)){
	nd <- append(nd, NA, which(colnames(d)==factors[i])-1)}
	names(nd)[which(colnames(d) %in% factors)]<-factors}

# if (!is.null(other) & any(!(names(other) %in% coefs))) {
	# wrong <- names(other)[which(!(names(other) %in% coefs))]
	# ifelse (length(wrong)>1, stop (paste("Variables", wrong, "not found in model")),
							 # stop (paste("Variable", wrong, "not found in model"))
			# )              
# }
if (!is.null(other)) {
	ph <- other[which(names(other) %in% coefs)]
	nd [match(names(ph), coefs)] <- as.vector(unlist(ph))
}

if(!is.null(values)) n <- length(values)

nd <- data.frame(matrix(nd,n,length(nd),byrow=T))
nd[,which(!(colnames(d) %in% factors))] <- unlist(lapply(as.list(nd)[which(!(colnames(d) %in% factors))],function(.)as.numeric(as.character(.))))
nd[,which(colnames(d) %in% factors)] <- unlist(lapply(data.frame(nd[which(colnames(d) %in% factors)]),as.character))
colnames (nd) <- coefs

xvalues <-seq(range[1],range[2],length.out=n)
if(!is.null(values)) xvalues <- values
nd [,grep(x,coefs)] <- xvalues

# This part is just so model.matrix wont get freaked out about full rank
if(length(factors)>0){
	extras <- prod(unlist(lapply(as.list(d)[factors],function(.)length(unique(.)))))
	nd <- rbind (nd, nd[(n-extras+1):n,])
	lvls <- lapply(as.list(d)[factors],function(.)levels(as.factor(.)))
	if (length(lvls)>1){
		for(i in 1:length(lvls)) nd[(n+1):(n + extras),names(lvls)[i]] <-as.vector(unlist(lvls[i]))
		} else nd[(n+1):(n + extras),names(lvls)] <- as.vector(unlist(lvls))
}

	
	if(use.predict) {
if(class(object)[1] %in% c("lmerMod","glmmPQL","glmerMod","lmerModLmerTest")){
ranints <- names(other)[which(names(other) %in% names(ranef(object)))]
if(length(ranints>0)) nd[,ranints] <- as.data.frame(lapply(other[ranints], function(x)rep(x, nrow(nd))))
}
nd <- nd[1:n,]
pr <- predict(object, newdata=nd)
return(pr)
}else{
nd <-   as.list(nd)

mm <- model.matrix(as.formula(paste("~",paste(vars, collapse="+"))), data=nd)
mm <- mm[1:n,]
if(class(object)[1] %in% c("lmerMod","glmmPQL","glmerMod","lmerModLmerTest")){
ranints <- names(other)[which(names(other) %in% names(ranef(object)))]
if(length(ranints>0)){
   for (i in 1:length(ranints)){
        opts <- ranef(object)[[match(ranints[i],names(ranef(object)))]]
		multiplier <- unlist(opts)[match(other[ranints][[i]], rownames(opts))]
 mm[,1] <- mm[,1]*multiplier
}}}
nd <- lapply(nd, function(x) x[1:n])

if (mixture){
  mmcond <- mm[,names(betas$cond)]
  mmzi <- mm[,names(betas$zi)]
  y_hat<-exp( mmcond%*%betas$cond)*(exp(mmzi%*%-betas$zi)/(1+exp(mmzi%*%-betas$zi)))
}else {y_hat<-mm%*%betas}

pr <- data.frame(cbind(xvalues,y_hat))
colnames (pr) [2] <- "predicted.mean"

if (CI) {
if(class(object)[1]=="lmerMod" | class(object)[1]=="glmerMod"){ predFun<-function(.) mm%*%fixef(.)
bb<-bootMer(object,FUN=predFun,nsim=CIsims,.progress=progress)
bb_se<-apply(bb$t,2,function(x) x[order(x)])
lo <- which.min(abs(1:CIsims - (100-conf)*CIsims/200))
hi <- CIsims - lo
pr$LowerCI <- bb_se[lo,]
pr$UpperCI <- bb_se[hi,]
}

if(class(object)[1]=="lm" | class(object)[1]=="glm"){
v <-vcov(object)
var.pred <- rowSums((mm %*% v) * mm) # more efficient way of calculating diag(mm %*% v %*% t(mm))
# this is because var(y_hat) = var(X*B_hat) = X*var(B_hat)*t(X)
se.pred<-sqrt(var.pred)
tval <- qt((100-conf)/200, df=(df.residual(object)))
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
  if(class(object)[1]=="glmmTMB"){
    if(mixture){
      params <- c(names(betas$cond),paste0("zi~",names(betas$zi)))
      vc <- vcov(object,full=T)[params,params]
      sim <- matrix(NA, CIsims, n)
      for (i in 1: CIsims){
        betas.sim <- rmvnorm(n=1,mean=c(betas$cond,betas$zi),sigma=vc)
        sim[i,] <- exp( mmcond%*%betas.sim[1:length(betas$cond)])*(exp(mmzi%*%-betas.sim[(length(betas$cond)+1):length(betas.sim)])/(1+exp(mmzi%*%-betas.sim[(length(betas$cond)+1):length(betas.sim)])))
      }
      pr$LowerCI <- apply(sim,2,quantile,(1-conf/100)/2)
      pr$UpperCI <- apply(sim,2,quantile,1-((1-conf/100)/2))
      
    }
  }
if(class(object)[1]=="brmsfit"){
    probs <- c((1-conf/100)/2,1-((1-conf/100)/2))
    brms_predict <- as.data.frame(predict(object, newdata=nd, probs=probs, re_formula=NA))
	pr$LowerCI <- brms_predict[,paste0("Q",probs[1]*100)]
    pr$UpperCI <- brms_predict[,paste0("Q",probs[2]*100)]
	pr$brms_estimate_check <- brms_predict$Estimate
}


}

if(unlink==T & !mixture){
pr[,2:ncol(pr)] <- family(object)$linkinv(as.matrix(pr[,2:ncol(pr)]))
}



return (pr)

}
}
