index <- function(data,row.index,col.index,na.rm=F){
args <- as.list(match.call())
#print(args)}
ph <- 'c('
if(na.rm) ph <- 'which('
data.name <- as.character(args[which(names(args)=='data')])
if('row.index' %in% names(args)){
row.index <- eval(parse(text=paste0('with(',data.name,',',ph,noquote(as.character(args[which(names(args)=='row.index')])),'))')))
}
if('col.index' %in% names(args)){
col.index <- eval(parse(text=noquote(as.character(args[which(names(args)=='col.index')]))))

}

data <- data[row.index,col.index]
return(data)
}