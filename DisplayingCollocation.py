from DESolver5 import *

mpmath.mp.dps = 100
mpmath.mp.prec = 100


# Defining the Domain and Number of Collocation Points
ps = 100
D = [0,1]
a = D[0]
b = D[1]
dx = (b-a)/(ps-1)
xlist = [a+i*dx for i in range(ps)]

# Defining our ODE
def f(x):
	return 0

def g1(x):

	return 1

def g0(x):

	return 1

g = [g0, g1]
XY = [[0],[1]]
N = 10

p = DESolver(g,f,D=D,XY=XY,DXDY=[],basis="c",N=N)

point = points(N-1,a,b)

def poly(x):
	return polyeval(p,x)

ay = [poly(x) for x in xlist]
ey = [math.e**(-x) for x in xlist]

plt.clf()
plt.title("Plot of $y=y_{10}$ versus $y=e^{-x}$")
plt.grid()
plt.plot(xlist,ey, "-b",label="$y=e^{-x}$")
plt.plot(xlist,ay, "--g",label="$y=y_{10}$")

plt.plot(0,1,"or",label="(0,1)")
plt.legend()
plt.show()

plt.clf()
def diffeq(x):
	return polyeval(polydiff(p.copy()),x)+polyeval(p,x)

ay = [diffeq(x) for x in xlist]
ey = [0 for x in xlist]

plt.grid()
plt.title("Putting $y_{10}$ into the ODE then plotting the LHS vs RHS")
plt.plot(xlist,ay, "-g",label="$y=LHS$")
plt.plot(xlist,ey, "-b",label="$y=RHS$")
i = 0
for po in point:
	if i == 0:
		plt.plot(po,0,"or",label="Collocation Points")
	else:
		plt.plot(po,0,"or")
	i += 1

plt.legend()
plt.show()

