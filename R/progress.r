progress <- function (i,length, percent=1) {
	if(i==1){print(0)}
	if(trunc(i*100/length) %% percent ==0 & trunc(i*100/length)!=trunc((i-1)*100/length)){ print(trunc(i*100/length))}
}
