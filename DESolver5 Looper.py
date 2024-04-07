from DESolver5 import *
from datetime import datetime


mpmath.mp.dps = 200
mpmath.mp.prec = 200

#  This procedure goes as follows:
#  We use Chebyshev to find a good solution on [0,1]
#  Then we sub x=1/u^k into our ODE
#  This transforms the [1,infinity] domain to [1,0] and we use Chebyshev to solve for g(u)=y(1/u^k) then using u=1/x we find y(x)=g(x^(-1/k))

# NOTE: For finding derivatives we can't use polydiff(g,x) for y(x) as y(x) is not a polynomial anymore...
plt.clf()
plt.grid()
plt.title("Pray This Works")
N = 25
it1 = 3
it2 = 6
n = 5

# Desired Domain
D = [0,25]


# First we use Chebyshev to create a good approximation on [0,1]

# Domain and Initial Guess
D1 = [0,1]
z1 = [1]

### ODE Solving
start=datetime.now()

# Iterative loop carrying out the QLM
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
	z1 = DESolver(g,f,D=D1,XY=[[0],[1]],DXDY=[[0],[0]],basis="c",N=N)


# Chooses the power for our sub
k = 1

# Defines the interval that we are finding a solution on
D2 = [D[1]**(-1/k),1]

# Creates our Initial Conditions
a2 = polyeval(z1,1)
b2 = -k*polyeval(polydiff(z1.copy()),1)

# Initial guess that satisfies our boundary conditions as well as decomposes as x goes to 0 
# (which reflects y goes to 0 as x goes to infinity before the substitution)
z2 = [0,0,3*a2-b2,b2-2*a2]

# Iterative loop carrying out the QLM
for i in range(it2):

	# OUr ODE's Coefficients
	def g0(x):
		return n*polyeval(z2,x)**(n-1)

	def g1(x):
		return (1-k)/k**2*x**(2*k+1)

	def g2(x):
		return x**(2*k+2)/k**2

	# What our ODE equals
	def f(x):
		return (n-1)*polyeval(z2,x)**n

	g = [g0,g1,g2]
	z2 = DESolver(g,f,D=D2,XY=[[1],[a2]],DXDY=[[1],[b2]],basis="c",N=N)

def f(x):
	if x < 1:
		return polyeval(z1, x)
	else:
		return polyeval(z2, x**(-1/k))

def g(x):
	return polyeval(z1,x)

points = 300
dx = (D[1]-D[0])/(points-1)
xlist = [D[0]+i*dx for i in range(points)]

print("Solution Found in: ", datetime.now()-start)

### Error Calculation
start = datetime.now()


# Calculates Residual Error for the solution: z1
def er1(x):
	return abs(x*polyeval(polydiffn(z1.copy(),2),x)+2*polyeval(polydiff(z1.copy()),x)+x*polyeval(z1.copy(),x)**n)

# Calculates Residual Error for the solution: z2
def er2(x):
	u = x**(-1/k)
	return abs(1/(k**2)*x**(-2/k-1)*polyeval(polydiffn(z2.copy(),2),u)+(1/k**2-1/k)*x**(-1/k-1)*polyeval(polydiff(z2.copy()),u)+x*(polyeval(z2,u)**n))

def er(x):
	if x < 1:
		return er1(x)
	else:
		return er2(x)

Error = integrate(er,D[0],D[1])

print("Raw Error: ", Error)
print("Error Normalised and Log'ed", -math.log(1/(D[1]-D[0])*Error,10))


print("Error Calculations Found in: ", datetime.now()-start)

### Plotting
start=datetime.now()
ylist = [f(x) for x in xlist]
y1list = [g(x) for x in xlist]
plt.plot(xlist,ylist, "-r",label="Final Solution")
ax = plt.gca()

plt.xlabel("x")
plt.ylabel("y")
plt.title("Lane-Emden Solution $n=3$")
plt.legend(loc="upper right")
ax.set_ylim([-0.3, 1.1])
print("Plots Created in: ", datetime.now()-start)

plt.show()


