from DESolver5 import *

mpmath.mp.dps = 200
mpmath.mp.prec = 200

M=5
n = M


plt.title(r"Chebyshev and rational Chebyshev domain patching (Lane-Emden $n=5$ on $[0,50]$)")

# Total Domain [0,L1]
D = [0,50]
L = 10

# For plots
points = 100
dx = (D[1]-D[0])/(points-1)
xlist = [D[0]+i*dx for i in range(points)]
plt.plot(xlist,[1/(sqrt(1+x**2/3)) for x in xlist],"-r",label=r"$y=\frac{1}{\sqrt{1+x^2/3}}$")

# Domain Chebyshev
D1=[0,L]

# Initial Guess Chebyshev
z = [1]
for i in range(10):

	def g0(x):
		return M*polyeval(z,x)**(M-1)

	def g1(x):
		return 2/x

	def g2(x):
		return 1

	# What our ODE equals
	def f(x):

		return (M-1)*polyeval(z,x)**(M)

	g = [g0,g1,g2]
	yi1 = z.copy()
	
	z = DESolver(g,f,D=D1,XY=[[0],[1]],DXDY=[[0],[0]],basis="c",N=25)


# Creating initial conditions for the next patch
b1 = polyeval(z,10)
b2 = polyeval(polydiff(z),10)
N1 = 25
it1 = 8

# Domain for Rat Cheby in x
D3 = [L,D[1]]

# Domain for Rat Cheby in u
D2 = [1/(1+D[1]),1/(1+D3[0])]


# Fix initial guess
z1 = [0,(1+L)*(2*b1+(1+L)*b2),-(1+L)**2*((1+L)*b2+b1)]

for i in range(it1):

	# Our ODE's Coefficients
	def g0(x):
		return n*polyeval(z1,x)**(n-1)

	def g1(x):
		return 2*x**4/(x-1)

	def g2(x):
		return x**4

	# What our ODE equals
	def f(x):
		return (n-1)*polyeval(z1,x)**n

	g = [g0,g1,g2]
	z2 = DESolver(g,f,D=D,XY=[[D2[1]],[b1]],DXDY=[[D2[1]],[-(1+L)**2*b2]],basis="c",N=N1)

	z1 = z2.copy()

def h(x):
	if x < L:
		return polyeval(z,x)
	else:
		return polyeval(z1,1/(1+x))

ylist = [h(x) for x in xlist]
plt.plot(xlist,ylist,"-g",label=r"$y=y_{approx}$")
ax = plt.gca()
plt.xlabel("x")
plt.ylabel("y")
plt.legend()

ax.set_ylim([0, 1.1])
plt.show()
plt.clf()