#library("kdensity")
library(cubature)
kdensity = function(x, bw) 
{
	f= function(z)
	{
		n=length(x)
		res =0
		for(i in 1:n)
		{
			res = res +K_g((z-x[i])/bw)
		}
		return(res)
	}
	return(f)
}
K_g = function(x)
{
	return(dnorm(x))
}
K_star = function(x)
{
	res = dnorm(x,0,sqrt(2)) - 2*K_g(x)
}
J = function(h,x)
{
	f = kdensity(x,bw =h)
	f2 = function(x) return(f(x)^2)
	res = cubintegrate(f2,lower = -5, upper = 5)$integral
	n = length(x)
	for(i in 1:n)
	{
		f=kdensity(x[-i],bw=h)
		res = res - 2*f(x[i])/n
	}
	return(res)
}

J_tilde = function(h,x)
{
	res=0
	n = length(x)
	for(i in 1:n)
	{
		for(j in 1:n)
		{
			res = res + K_star((x[i]-x[j])/h)
		}
	}
	res = res/(n^2*h)
	res = res + 2*K_g(0)/(n*h)
	return(res)
}

J_star = function(h,x)
{
	f = kdensity(x,bw =h)
	f2 = function(x) return(f(x)^2)
	res = cubintegrate(f2,lower = -Inf, upper = Inf)$integral
	n = length(x)
	for(i in 1:n)
	{
		res = res - 2*f(x[i])/n
	}		
	return(res)
}

Run_J = function(x,mu,sd)
{
	fun = function(h)
	{
		return(J(h,x))
	}
	#print(fun(0.5))
	n=length(x)
	h_J = optim(0.5,fun,method ="Brent",lower =0 , upper =10*n^(-1/5))$par
	#print(h_J)
	z= density(x,bw = h_J , kernel ="gaussian")
	df = approxfun(z)
	f = function(x)
	{
		return((df(x)-dnorm(x,mu,sd))^2)
	}
	#print(f(1))
	ise = cubintegrate(f,lower = max(-5,min(z$x)), upper = min(5,max(z$x)))$integral
	#print(ise)
	return(ise)
}

Imse_J = function(iter=1000,n=100,mu =0 , sd =1)
{
	x = rep(n,iter)
	#print(x[1])
	f = function(x) return(rnorm(x,mean = mu,sd = sd))
	z= lapply(as.list(x),f)
	#print(z[[1]])
	f = function(x)
	{
		return(Run_J(x,mu,sd))
	}
	res=lapply(z,f)
	#print(res)
	return(Reduce('+',res)/iter)

}

Run_J_tilde = function(x,mu,sd)
{
	fun = function(h)
	{
		return(J_tilde(h,x))
	}
	#print(fun(0.5))
	n=length(x)
	h_J = optim(0.5,fun,method ="Brent",lower =0 , upper =10*n^(-1/5))$par
	#print(h_J)
	z= density(x,bw = h_J , kernel ="gaussian")
	df = approxfun(z)
	f = function(x)
	{
		return((df(x)-dnorm(x,mu,sd))^2)
	}
	#print(f(1))
	ise = cubintegrate(f,lower = max(-5,min(z$x)), upper = min(5,max(z$x)))$integral
	#print(ise)
	return(ise)
}
Imse_J_tilde = function(iter=1000,n=100,mu =0 , sd =1)
{
	x = rep(n,iter)
	#print(x[1])
	f = function(x) return(rnorm(x,mean = mu,sd = sd))
	z= lapply(as.list(x),f)
	#print(z[[1]])
	f = function(x)
	{
		return(Run_J_tilde(x,mu,sd))
	}
	res=lapply(z,f)
	#print(res)
	return(Reduce('+',res)/iter)
}
Run_J_star = function(x,mu,sd)
{
	fun = function(h)
	{
		return(J_star(h,x))
	}
	#print(fun(0.5))
	n=length(x)
	h_J = optim(0.5,fun,method ="Brent",lower =0 , upper =10*n^(-1/5))$par
	#print(h_J)
	z= density(x,bw = h_J , kernel ="gaussian")
	df = approxfun(z)
	f = function(x)
	{
		return((df(x)-dnorm(x,mu,sd))^2)
	}
	#print(f(1))
	ise = cubintegrate(f,lower = max(-5,min(z$x)), upper = min(5,max(z$x)))$integral
	#print(ise)
	return(ise)
}
Imse_J_star = function(iter=1000,n=100,mu =0 , sd =1)
{
	x = rep(n,iter)
	#print(x[1])
	f = function(x) return(rnorm(x,mean = mu,sd = sd))
	z= lapply(as.list(x),f)
	#print(z[[1]])
	f = function(x)
	{
		return(Run_J_star(x,mu,sd))
	}
	res=lapply(z,f)
	#print(res)
	return(Reduce('+',res)/iter)
}
print(Imse_J_star(iter =30))