from DESolver5 import *

mpmath.mp.dps = 500
mpmath.mp.prec = 500


plt.clf()

plt.grid()
plt.xlabel(r"$N$ Value")
plt.ylabel("Residual Error")
plt.title(r"Plot of Residual Error for the Lane Emden Equation with $n=3$")
K = 51
D = [0,10]
xlist = []
ylist = []

def g0(x):
		return 3

def g1(x):
	return 2/x

def g2(x):
	return 1

# What our ODE equals
def f(x):
	return 2

g = [g0,g1,g2]

# OUr ODE's Coefficients
for N in range(3,K,2):
	
	z2 = DESolver(g,f,D=D,XY=[[0],[1]],DXDY=[[0],[0]],basis="c",N=N)

	def ser(x):
		return abs(x*3*polyeval(z2,x)+2*polyeval(polydiff(z2),x)+x*polyeval(polydiffn(z2,2),x)-2*x)

	E = -math.log(integrate(ser,0,D[1]),10)/D[1]
	print(E, "Residual Error ", N)
	print(10**(E*(-D[1])))
	xlist.append(N)
	ylist.append(E)

y1list = []
for N in range(3,K,2):

	z2 = DESolver(g,f,D=D,XY=[[0],[1]],DXDY=[[0],[0]],basis="t",N=N)
	def ser(x):
		return abs(x*3*trigeval(z2,x)+2*trigeval(trigdiff(z2),x)+x*trigeval(trigdiffn(z2,2),x)-2*x)

	E = -math.log(integrate(ser,0,D[1]),10)/D[1]
	print(E, "Residual Error ", N)
	print(10**(E*(-D[1])))

	y1list.append(E)


plt.plot(xlist,ylist,":r",label="Residual Error: Chebyshev Basis")
plt.plot(xlist,y1list,":g",label="Residual Error: Trigonometric Basis")

plt.legend()
plt.show()
