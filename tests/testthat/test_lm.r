

test_that("basic linear model works", {
	N<-25
	a<- 6
	x1<-runif(N,-100,100)
	b1<- 2
	x2<-rbinom(N,20,0.2)
	b2<-9
	x3<-rnorm(N,20,50)
	b3<- -4
	se<-50
	y<- a + b1*x1 + b2*x2 + b3*x3 +rnorm(N,0,se)
	d<- as.data.frame(cbind(y,x1,x2,x3))
	mod<-lm(y~x1+x2+x3,data=d)
	expect_silent(out <- counterfact(mod, x="x1", CI=T, other=list(x2=21)))
	expect_equal(nrow(out), 1000)
	expect_error(counterfact(mod, x="x1", CI=T, other=list(nonsense=21)))
	expect_equal(ncol(out),4)
	
})