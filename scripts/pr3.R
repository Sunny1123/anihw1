Run = function(n=200,h = 5,mu=0,sd =1,k = "gaussian")
{
	x= rnorm(n,mu,sd)
	z= density(x,bw = h , kernel =k)
	df = approxfun(z)
	f = function(x)
	{
		return((df(x)-dnorm(x,mu,sd))^2)
	}
	ise = integrate(f,lower = min(z$x), upper = max(z$x))$value
	return(ise)
}
Imse = function(iter=1000,n=200,mu =0 , sd =1,h = 5,k = "gaussian")
{
	x = rep(n,iter)
	f = function(x)
	{
		return(Run(x,mu =mu ,sd =  sd,h= h,k = k))
	}
	z=sapply(x,f)
	return(Reduce('+',z)/iter)

}

main = function(...)
{
	n= 10:100
	eff1= 0.9511
	eff2 = 0.9295
	f = function(x,...)
	{
		return(Imse(n=x,...))
	}
	epanechnikov = sapply(n,function(x) f(x,k = "epanechnikov"))
	gauss = sapply(floor(n/eff1),function(x) f(x,k = "gaussian"))
	boxcar = sapply(floor(n/eff2),function(x) f(x,k = "rectangular"))
	jpeg("pr3.jpg" , width = 1080 , height = 360)
	par(mfrow =c(1,3))
	plot(n ,epanechnikov , ty ='l', xlab = "Sample size")
	plot(floor(n/eff1), gauss, ty ='l', xlab = "Sample size")
	plot(floor(n/eff2), boxcar, ty = 'l',xlab = "Sample size")
	dev.off()
}
main(iter =100)
