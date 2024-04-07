from DESolver5 import *

v = 0.9
D = [-v,v]
def g(x):
	return 0

def f2(x):
	return 1-x**2

def f1(x):
	return -x

def f0(x):
	return 1/25

plt.title(r"$(1-x^2)y''-xy'+\frac{1}{25}y=0$: Approximate solutions $y_N$ for $N\in\{2,3,...,10\}$ on $D=[-"+str(v)+","+str(v)+"]$")
points = 100
dx = (D[1]-D[0])/points
xlist = [D[0]+i*dx for i in range(points+1)]
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")

f = [f0,f1,f2]
for i in range(3,27):

	p = DESolver(f,g,D,XY=[[0],[1]],DXDY=[[0],[1]], basis="c", N=i)

	ylist = [polyeval(p,x) for x in xlist]
	plt.plot(xlist,ylist,"--",label=r"$y_{"+str(i-1)+"}$")


	def err(x):
		return abs((1-x**2)*polyeval(polydiff(polydiff(p)),x)-x*polyeval(polydiff(p),x)+1/25*polyeval(p,x))

	print(-math.log(1/(2*v)*integrate(err,-v,v),10),i-1)


plt.legend()
plt.show()