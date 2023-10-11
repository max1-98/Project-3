

# In this again we will solve an N dimensional linear DE. 


from main import *
import math
import numpy as np
import matplotlib.pyplot as plt

### Inputs from the User

# Accuracy (number of Orthogonal Polys)
N = 3

# Coefficients of the terms in the ODE "Coeffunctions"
# From y, dy/dx, d^2y/dx^2...
def g0(x):
	return 1

def g1(x):
	return 1 

def g2(x):
	return 0

def g3(x):

	return x**4


# What our ODE equals
def f(x):

	return 0

# Initial Conditions

X = [0]
Y = [1]

# The solutions domain [a,b]
D = [-1,1]


# Our plan:

# First we will create and transform our Chebyshev nodes to [a,b]
# Then we will transform our Chebyshev Polynomials to [a,b]

# Then we will input our initial points to create n equations
# Finally we will do collocation on N-n points to create the remaining N-n equations

# Solve this via GJ

# Our coefficients
coeffun = [g0,g1]

# Order of the ODE
n = len(coeffun)-1


# Collocation Points
points1 = ChebyshevNodes(N-n)
a = D[0]
b = D[1]

# Transformed Points
points2 = []

for point in points1:
	points2.append((point*(b-a)/2+(b+a)/2))

# Original Chebyshev Polynomials
polysNT = ChebyshevGen(N)

#Creating Our Translated Versions
polys = []
for i in range(N):
	p1 = polysNT[i]
	p1 = polyxstret(p1,(b-a)/2)
	p1 = polyxtrans(p1,(b+a)/2)
	
	polys.append(p1)

# So polys are our transformed Chebyshev Polynomials

### Next let's create the N equations


# Mx=B where x is the vector representing our coefficients a0,a1,a2...
M = []
B = []

# We'll start with the Initial Conditions
for i in range(n):

	R=[]

	for j in range(N):

		R.append(polyeval(polys[j],X[i]))
	B.append(Y[i])
	M.append(R)


# Next we need our N-n equations from plugging in yN and evaluating it at our transformed Chebyshev Nodes

for i in range(N-n):

	# ith Transformed Chebyshev Node
	x = points2[i]

	# Initialising ith row of M
	R = []

	# Creating R
	for j in range(N):

		# Initialising the value of the coefficient of aj in the ith row
		S = 0

		# Copying our polynomial so it doesn't get editted by our functions
		p = polys[j].copy()


		# Calculating the coefficient of aj in the ith row
		for k in range(n+1):
			p1 = p.copy()
			p2 = polydiffn(p1,k)
			

			S += coeffun[k](x)*polyeval(p2,x)

		R.append(S)

	M.append(R)
	B.append(f(x))


M = np.array(M)
B = np.array(B)

a = np.linalg.solve(M,B)


solution = []
for i in range(N):
	solution = polyadd(solution,polysmult(polys[i],a[i]))


def ef(x):

	return math.e**(-x)


points = 100
dx = (b-a)/(points-1)
xlist = [a+i*dx for i in range(points)]
eylist = [ef(x) for x in xlist]
aylist = [polyeval(solution,x) for x in xlist]



plt.plot(xlist,eylist,"-r")
plt.plot(xlist,aylist,"-g")
plt.show()