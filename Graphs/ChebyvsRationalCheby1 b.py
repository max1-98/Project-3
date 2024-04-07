from DESolver5 import *

v = 100
D = [1/(v+1),1]
def g(x):
	return 0

def f1(x):
	return x**2

def f0(x):
	return -1

plt.title(r"Rational Chebyshev Function Solution $y_{25}$ to: $\frac{dy}{dx}=-y\quad y(0)=1\quad x\in[0,100]$")

plt.xlabel(r"$x$")
plt.ylabel(r"$y$")

f = [f0,f1]


p = DESolver(f,g,D,XY=[[1],[1]],DXDY=[], basis="c", N=25)

def h(x):
	return polyeval(p,1/(x+1))

D1= [0,v]
points = 300
dx = (D1[1]-D1[0])/points
xlist = [D1[0]+i*dx for i in range(points+1)]
aylist = [h(x) for x in xlist]
eylist = [math.e**(-x) for x in xlist]

def er(x):
	return abs(h(x)-1/(x+1)**2*polyeval(polydiff(p),1/(x+1)))

print(-math.log(1/100*integrate(er,0,100),10))


plt.plot(xlist,eylist,"-r",label=r"$y^*$")
plt.plot(xlist,aylist,"--b",label=r"$y_{25}$")



plt.legend()
plt.show()


















