from DESolver5 import *
from datetime import datetime

mpmath.mp.dps = 200
mpmath.mp.prec = 200


N = 25
it = 10
plt.clf()
plt.title(r"Using Quasilinearization to solve: $y(x)y'(x)=x\quad y(0)=1$ on $[-1,1]$")
plt.grid()

D = [-1,1]

points = 300
dx = (D[1]-D[0])/points
xlist = [dx*i+D[0] for i in range(points+1)]
eylist = [math.sqrt(x**2+1) for x in xlist]
plt.plot(xlist,eylist,"-b",label="Exact Solution")

z = [1]
sylist = [polyeval(z,x) for x in xlist]
plt.plot(xlist,sylist,"--",label="Iterate: 0")

start = datetime.now()

for i in range(it):
	def g0(x):
		return polyeval(polydiff(z),x)

	def g1(x):
		return polyeval(z,x)

	def f(x):
		return x+polyeval(polydiff(z),x)*polyeval(z,x)

	g = [g0,g1]

	z = DESolver(g,f,D=D,XY=[[0],[1]],DXDY=[],basis="c",N=N)

	ylist = [polyeval(z,x) for x in xlist]
	plt.plot(xlist,ylist,"--",label="Iterate: "+str(i+1))
	def h(x):
		return abs(polyeval(z,x)*polyeval(polydiff(z),x)-x)
	print(-math.log(1/(D[1]-D[0])*integrate(h,D[0],D[1]),10), " Iterate: "+str(i+1))

ax = plt.gca()
plt.xlabel("x")
plt.ylabel("y")
#ax.set_ylim([0, 100])
plt.legend()
plt.show()