progress <- function (i,length, increment=100,step=1) {
	if(i==1){print(0)}
	if(trunc(i*increment/length) %% step ==0 & trunc(i*increment/length)!=trunc((i-1)*increment/length)){ print(trunc(i*increment/length))}
if(i==length) print("Done")
}
