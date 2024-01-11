# This DESolver is not an advancement. This is for using Cardinal functions. Later versions we will create a DEsolver that allows us
# to use both. 

# This DE Solver has had the following modified as the basis of Cardinal Functions Naturally satisfied one of the initial conditions
# - N-n was changed to N-n+1
# - The first point was changed to x=1 instead of x=0 as this matched the initial value

### Required Imports
from main import *
from mpmath import *
import cmath
import math
import matplotlib.pyplot as plt



# Cardinal Function r with shifting constant a and width w
def f(x,r,a=0,w=1):

	if x == a+r*w:
		return 1
	else:
		return sin(pi/w*(x-a-r*w))/(pi/w *(x-a-r*w))

# The sum of all the cardinal functions with their respective coefficients
def G(x,A,a=0,w=1):

	s = 0
	for i in range(len(A)):
		s += A[i]*f(x,i,a,w)

	return s


### Defining our Linear n Order ODE

# coefficient of y

def g0(x):
	return 1
# coefficient of dy/dx
def g1(x):
	return 1

# coefficient of dy2/dx2
def g2(x):
	return x

# What our ODE equals
def h(x):

	return 0
"""
def g0(x):
	return x*7*G(x,A,a,w)**6

# coefficient of dy/dx
def g1(x):
	return 2

# coefficient of dy2/dx2
def g2(x):
	return x

# What our ODE equals
def h(x):

	return x*6*G(x,A,a,w)**7
"""

# Our Domain
D = [0,3]

# List of the functions
#g = [g0,g1,g2]
g = [g0,g1]

# Calculating the order of our ODE
n = len(g)-1

# Initial Conditions
"""
XY = [[0,1]]
DXDY = [[0,0]]
"""
XY = [[0,1]]
DXDY = []

# Selecting how many Cardinal Functions we want
N = 50

# Defining a and w
a = D[0]
w = (D[1]-D[0])/(N-n+1)



# Our Initial Guess
A = [0]*N
A[0] = 1


### --- Calculations --- ###

### Setup



# Now we know the order we can generate our Collocation points, calculate the derivatives ext.
points = []
for i in range(N-n+1):
	points.append(0.1+a+i*w)

# Here we create a triple-nested loop containing the values of our function and it's derivatives evaluated at the points found above

# In total (n+1)(N-n)(N) points in our list so our complexity is Quadratic -- Not Ideal but not too bad
KV = []
for i in range(n+1):
	# Create our Nested List
	S = []
	# Loops through the Collocation Points
	for j in range(N-n+1):
		# Loops through the Different Cardinal Functions

		R = []
		# Creates a row where we evaluate the kth Cardinal Function's ith derivative at the jth node
		for k in range(N):

			if i == 0:
				R.append(f(points[j],k,a,w))
			else:
				R.append(CSD(k,points[j],a,1,i,10**(-8)))

		S.append(R)
	KV.append(S)

for i in range(10):
	### Forming the Equations MA = B

	M = []
	B = []
	# Initial Conditions

	# Coordinates
	for i in range(len(XY)):

		B.append(XY[i][1])
		R = []
		for j in range(N):
			R.append(f(XY[i][0],j,a,w))

		M.append(R)

	# Derivatives

	"""
	for i in range(len(DXDY)):

		B.append(DXDY[i][1])
		R = []

		for j in range(N):
			R.append(CSD(j,DXDY[i][0],a,1,1,10**(-8)))
		print(R)
		M.append(R)
	"""

	# Equations from Collocation Points

	for i in range(N-n+1):
		R = []
		x = points[i]

		for j in range(N):
			s = 0
			for k in range(n+1):
				s += g[k](x)*KV[k][i][j]
			R.append(s)
		M.append(R)
		B.append(h(x))

	A = mpmath.lu_solve(M,B)

# Creating the plot
ps = 200
dx = (D[1]-D[0])/ps
xlist = [D[0]+dx*i for i in range(ps)]
ylist = [G(x,A,a,w) for x in xlist]
plt.plot(xlist,ylist, "-r")
plt.show()



