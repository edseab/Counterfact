unique_permuts <- function(vector, mirrors=T){
Var1 <- c()
Var2 <- c()

if(mirrors){
for(i in 1:(length(vector))){
Var1 <- c(Var1, rep(vector[i], (length(vector)-i+1)))
Var2 <- c(Var2, vector[(i):length(vector)])
}}

if(!mirrors){
for(i in 1:(length(vector)-1)){
Var1 <- c(Var1, rep(vector[i], (length(vector)-i)))
Var2 <- c(Var2, vector[(i+1):length(vector)])
}}
return(data.frame(Var1=Var1,Var2=Var2))
}