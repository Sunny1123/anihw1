library("kdensity")

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
	f = kdensity(x,bw =h,kernel ="gaussian")
	f2 = function(x) return(f(x)^2)
	res = integrate(f2,lower = -Inf, upper = Inf)$value
	n = length(x)
	for(i in 1:n)
	{
		res = res - 2*kdensity(x[-i],bw=h,kernel = "gaussian")(x[i])/n
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
			res += K_star((x[i]-x[j])/h)
		}
	}
	res = res/(n^2*h)
	res += 2*K_g(0)/(n*h)
	return(res)
}
