from DESolver5 import *

def f(x):
	return 0


def g0(x):
	return 1/25

def g1(x):
	return (x+1)*(3*x+2)

def g2(x):
	return (x+1)**2*(1+2*x)

D1 = [0, 0.95]
D = [0,1+1/(1-D1[1])]
N=25




plt.title(r"Rational Chebyshev Function Solution $y_{25}$ to: $\frac{dy}{dx}=-y\quad y(0)=1\quad x\in[0,100]$")

plt.xlabel(r"$x$")
plt.ylabel(r"$y$")

g = [g0,g1,g2]


p = DESolver(g,f,D,XY=[[0],[1]],DXDY=[[0],[1]], basis="l", N=N)

def h(x):
	return polyeval(p,-1-1/(x-1))

points = 300
dx = (D1[1]-D1[0])/points
xlist = [D1[0]+i*dx for i in range(points+1)]
aylist = [h(x) for x in xlist]

plt.plot(xlist,aylist,"--b",label=r"$y_{25}$")

def f(x):
	return 0

def g0(x):
	return 1/25

def g1(x):
	return -1*x

def g2(x):
	return (1-x**2)

g = [g0,g1,g2]


p = DESolver(g,f,D1,XY=[[0],[1]],DXDY=[[0],[1]], basis="c", N=N)

plt.plot(xlist,[polyeval(p,x) for x in xlist],"--r")
plt.legend()
plt.show()