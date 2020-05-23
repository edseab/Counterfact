w2l <- function (data, vars, sep=""){
repeat_number <- ncol(data)
full_cols <- paste(rep(vars, each=repeat_number), rep(1:repeat_number, length(vars)), sep=sep)
options <- colnames(data)[colnames(data) %in% full_cols]
newcols <- matrix (NA,nrow=nrow(data), ncol=sum(!full_cols %in% colnames(data)))
colnames(newcols) <- full_cols[(!full_cols %in% colnames(data))]
data <- cbind(data, newcols)
new_data <- do.call("rbind", replicate(repeat_number, data[!colnames(data) %in% full_cols], simplify = FALSE))
for (i in 1:length(vars)){
new_var <- stack(data[,paste(vars[i],1:repeat_number,sep=sep)])
new_data <- cbind(new_data,new_var$values)
colnames(new_data)[ncol(new_data)] <- vars[i]
}
new_data <- new_data[which(!apply(new_data[,vars],1,function(x) all(is.na(x)))),]
return(new_data)
}
