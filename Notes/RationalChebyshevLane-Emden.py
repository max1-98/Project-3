from DESolver5 import *

mpmath.mp.dps = 100
mpmath.mp.prec = 100

#  This procedure goes as follows:
#  We use Chebyshev to find a good solution on [0,1]
#  Then we sub x=1/u into our ODE
#  This transforms the [1,infinity] domain to [1,0] and we use Chebyshev to solve for g(u)=y(1/u) then using u=1/x we find y(x)=g(1/x)

# NOTE: For finding derivatives we can't use polydiff(g,x) for y(x) as y(x) is not a polynomial anymore...

plt.clf()
plt.grid()
plt.title("Pray This Works")
N = 25
it1 = 5
it2 = 10
n = 3

# First we use Chebyshev to create a good approximation on [0,1]

D = [0,1]
z1 = [1]
for i in range(it1):

	# OUr ODE's Coefficients
	def g0(x):
		return n*polyeval(z1,x)**(n-1)

	def g1(x):
		return 2/x

	def g2(x):
		return 1

	# What our ODE equals
	def f(x):
		return (n-1)*polyeval(z1,x)**n

	g = [g0,g1,g2]
	z1 = DESolver(g,f,D=D,XY=[[0],[1]],DXDY=[[0],[0]],basis="c",N=N)

def er1(x):
	return abs(x*polyeval(polydiffn(z1.copy(),2),x)+2*polyeval(polydiff(z1.copy()),x)+polyeval(z1.copy(),x)**n)

print(-math.log(integrate(er1,0,D[1])/D[1],10))

a = polyeval(z1,1)
b = -polyeval(polydiff(z1.copy()),1)

D = [1/25,1]

z2 = [0,a-b+1,b-2,1]
for i in range(it2):

	# OUr ODE's Coefficients
	def g0(x):
		return n*polyeval(z2,x)**(n-1)

	def g1(x):
		return 0

	def g2(x):
		return x**4

	# What our ODE equals
	def f(x):
		return (n-1)*polyeval(z2,x)**n

	g = [g0,g1,g2]
	z2 = DESolver(g,f,D=D,XY=[[1],[a]],DXDY=[[1],[b]],basis="c",N=N)

def er2(x):
	return abs(x**4*polyeval(polydiffn(z2.copy(),2),x)+polyeval(z2,x)**n)

def er3(x):
	return abs(x*(1/x**4*polyeval(polydiffn(z2.copy(),2),1/x)+2/x**3*polyeval(polydiff(z2.copy()),1/x))+2/x**2*polyeval(polydiff(z2.copy()),1/x)+polyeval(z2,1/x)**n)

print(-math.log(integrate(er2,D[0],D[1])/D[0],10))
#print(-math.log(integrate(er3,1/D[1],1/D[0])/(-1/(D[1])+1/(D[0])),10))

def f(x):
	if x < 1:
		return polyeval(z1, x)
	else:
		return polyeval(z2, 1/x)

def g(x):
	return polyeval(z1,x)

print(z1,"z1")
print(z2,"z2")

D = [0,1/D[0]]
points = 200
dx = (D[1]-D[0])/(points-1)
xlist = [D[0]+i*dx for i in range(points)]

ylist = [f(x) for x in xlist]
y1list = [g(x) for x in xlist]
plt.plot(xlist,y1list,"--g")
plt.plot(xlist,ylist, "-r")
plt.plot(1,a,"oy")
#plt.legend()
ax = plt.gca()

ax.set_ylim([-1, 1.3])
plt.show()