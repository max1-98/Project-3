from DESolver5 import *

mpmath.mp.dps = 300
mpmath.mp.prec = 300

M=5
n = M


plt.title(r"Chebyshev and rational Chebyshev domain patching (Lane-Emden $n=5$ on $[0,50]$)")

# Total Domain [0,L1]
D = [0,50]

## Patches

# Chebyshev
c = 10

# Rational Chebyshev 
rc = 10

# For plots
points = 100
dx = (D[1]-D[0])/(points-1)
xlist = [D[0]+i*dx for i in range(points)]
plt.plot(xlist,[1/(sqrt(1+x**2/3)) for x in xlist],"-r",label=r"$y=\frac{1}{\sqrt{1+x^2/3}}$")

# Domain Chebyshev

D1=[0,10]

## Equally wide patches

pw = 10/c

# List of functions on each patch 
zp = []
# patches are [0,pw], [pw,2pw], [2pw,3pw]...[(c-1)*pw,c*pw]

# Initial guess Chebyshev patch 1
z = [1]
b1 = 1
b2 = 0 

for j in range(c):
	Dc = [j*pw,(j+1)*pw]
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
		
		z = DESolver(g,f,D=Dc,XY=[[Dc[0]],[b1]],DXDY=[[Dc[0]],[b2]],basis="c",N=25)

	# Creating initial conditions for the next patch
	b1 = polyeval(z,Dc[1])
	b2 = polyeval(polydiff(z),Dc[1])

	# Adds z to the list of z's 
	zp.append(z.copy())

	# Creates next initial guess
	z = [b1-b2*Dc[1],b2]

N1 = 25
it1 = 8

# Domain for Rat Cheby in x
D3 = [Dc[1],D[1]]

# Domain for Rat Cheby in u
D2 = [1/(1+D[1]),1/(1+D3[0])]

pw2 = (D2[1]-D2[0])/rc


# Initial guess first patch
z1 = [0,(1+Dc[1])*(2*b1+(1+Dc[1])*b2),-(1+Dc[1])**2*((1+Dc[1])*b2+b1)]

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
	z2 = DESolver(g,f,D=D,XY=[[D2[1]],[b1]],DXDY=[[D2[1]],[-(1+Dc[1])**2*b2]],basis="c",N=N1)

	z1 = z2.copy()

def h(x):
	if x < Dc[1]:
		i = int(x/pw)
		return polyeval(zp[i],x)
	else:
		return polyeval(z1,1/(1+x))

def err1(x):
	i = int(x/pw)

	if abs(i-x/pw)< 0.000000000001:
		i+=-1

	return abs(x*polyeval(polydiffn(zp[i],2),x)+2*polyeval(polydiff(zp[i]),x)+x*polyeval(zp[i],x)**M)

def err2(x):
	return abs(x**4*polyeval(polydiffn(z1,2),x)+2*x**4/(x-1)*polyeval(polydiff(z1),x)+polyeval(z1,x)**M)*1/(x)**2

r1 = integrate(err1,0,Dc[1])
r2 = integrate(err2,1/(1+D[1]),1/(1+Dc[1]))
print(-math.log(1/(D[1]-D[0])*(r1+r2),10))

ylist = [h(x) for x in xlist]
plt.plot(xlist,ylist,"-g",label=r"$y=y_{approx}$")
ax = plt.gca()
plt.xlabel("x")
plt.ylabel("y")
plt.legend()

ax.set_ylim([0, 1.1])
plt.show()
plt.clf()