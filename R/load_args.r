load_args <- function(ph,func=NULL){
ph[!sapply(ph,is.numeric)] <- paste0("'",ph[!sapply(ph,is.numeric)],"'")
if(!is.null(func))ph <- ph[names(ph) %in% names(formals(func))]
for(i in 1:length(ph)){eval(parse(text=paste(names(ph),ph,sep="<-")[i]),envir=.GlobalEnv)
print(paste(names(ph),ph,sep="<-")[i])}
print("Done")

}