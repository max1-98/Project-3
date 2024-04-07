from DESolver5 import *

D = [0,1]
def g(x):
	return 0

def f1(x):
	return 1

def f0(x):
	return -1

plt.title(r"Exact Solution of  $\frac{dy}{dx}=y$ versus approximate solutions $y_N$ for $N\in\{1,2,3\}$")
points = 100
dx = (D[1]-D[0])/points
xlist = [D[0]+i*dx for i in range(points+1)]
ylist = [math.e**x for x in xlist]
plt.plot(xlist,ylist,"-r",label=r"$y=e^x$")

f = [f0,f1]
for i in range(2,5):

	p = DESolver(f,g,D,XY=[[0],[1]],DXDY=[], basis="c", N=i)

	ylist = [polyeval(p,x) for x in xlist]
	plt.plot(xlist,ylist,"--",label=r"$y_{"+str(i-1)+"}$")

	def er(x):
		return abs(math.e**x-polyeval(p,x))
	print(-math.log(integrate(er,0,1),10),i-1)

	def err(x):
		return abs(polyeval(polydiff(p),x)-polyeval(p,x))

	print(-math.log(integrate(err,0,1),10),i-1)

plt.legend()
plt.show()