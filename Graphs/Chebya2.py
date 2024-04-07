from DESolver5 import *

D = [-1,1]
def g(x):
	return 0

def f2(x):
	return 1-x**2

def f1(x):
	return -x

def f0(x):
	return 4

plt.title(r"Exact Solution of  $(1-x^2)y''-xy'+4y=0$ versus approximate solutions $y_N$ for $N\in\{2,3,4\}$")
points = 100
dx = (D[1]-D[0])/points
xlist = [D[0]+i*dx for i in range(points+1)]
ylist = [1-2*x**2 for x in xlist]
plt.plot(xlist,ylist,"-r",label=r"$y=1-2x^2$")
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")

f = [f0,f1,f2]
for i in range(3,6):

	p = DESolver(f,g,D,XY=[[0],[1]],DXDY=[[0],[0]], basis="c", N=i)

	ylist = [polyeval(p,x) for x in xlist]
	plt.plot(xlist,ylist,"--",label=r"$y_{"+str(i-1)+"}$")

	"""
	def er(x):
		return abs(math.e**x-polyeval(p,x))
	print(-math.log(integrate(er,0,1),10),i-1)

	def err(x):
		return abs(polyeval(polydiff(p),x)-polyeval(p,x))

	print(-math.log(integrate(err,0,1),10),i-1)
	"""

plt.legend()
plt.show()
