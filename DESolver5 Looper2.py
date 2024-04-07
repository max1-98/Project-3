from DESolver5 import *

mpmath.mp.dps = 100
mpmath.mp.prec = 100

k = 1
D = [0,10]
N = 25
n = 5
it = 13
z = [0,0,3,-2]
D1 = [1/(D[1]+1)**k,1/(D[0]+1)**k]

for i in range(it):

	# OUr ODE's Coefficients
	def g0(x):
		return n*polyeval(z,x)**(n-1)

	def g1(x):
		return k*x**(2/k+1)*(1+k+2/(x**(1/k)-1))

	def g2(x):
		return k**2*x**(2/k+2)

	# What our ODE equals
	def f(x):
		return (n-1)*polyeval(z,x)**n

	g = [g0,g1,g2]
	z = DESolver(g,f,D=D1,XY=[[1],[1]],DXDY=[[1],[0]],basis="c",N=N)

def f(x):
	return polyeval(z,1/(x+1)**k)


points = 500
dx = (D[1]-D[0])/points
xlist = [D[0]+i*dx for i in range(points+1)]
ylist = [f(x) for x in xlist]
plt.plot(xlist,ylist,"-r")
plt.show()