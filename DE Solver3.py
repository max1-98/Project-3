
# This solver will solve any n order LINEAR differential equation


# Imports 
import numpy as np
import math
import matplotlib.pyplot as plt

# Inputs

# How accurate, for now accuracy will just be the number of orthogonal polynomials, later we will create more accuracy by adding greater precision
N = 10

# coefficients in the linear DE, NOTE: The 0th index represents coefficient of y, the nth index represents the nth derivative of y.
a = [0,-1,0,1]

# So our order is
o = len(a)-1

# Dy=f

# our function in the differential equation Dy=f where y is our solution and D is a linear differential operator. 

def f(x):
	return 0

# Boundary Values

X = [0,1,2]
Y = [-1,0,1]


# In this model we'll only do Chebyshev. So they will pick and then the following will be generated:

def ChebyshevGen(N):

	# We will create N Chebyshev polynomials

	C = []
	
	for i in range(2):
		T = [0]*N
		T[i] = 1
		C.append(T)

	for i in range(N-2):

		# First we shift the previous polynomial and multiply by 2.
		T = [0]*N
		Y = C[i+1]

		for j in range(N-1):
			T[j+1] = 2*Y[j]


		# Then we term wise subtract to the one before that. 
		Z = C[i]
		for j in range(N):
			T[j] = T[j] - Z[j]

		C.append(T)

	return C

# We need a way of differentiating 

def deriv(x):

	# x is a list of coefficient of the polynomial a0+a1x^1...
	l = len(x)
	y = [0]*len(x)

	for i in range(l-1):

		y[i] = (i+1)*x[i+1]

	return y

# and evaluating polynomials

def pe(x,a):
	# Polynomial Eval
	n = len(a)

	S = 0 
	for i in range(n):

		S += a[i]*x**i

	return S

# Finally we will calculate the derivatives of Chebyshevs that we need

"""(Chebyshev Derivatives)"""
TDs = [ChebyshevGen(N)]


for i in range(len(a)-1):
	TD = []
	TDI = TDs[i]
	for j in range(N):

		TD.append(deriv(TDI[j]))
	TDs.append(TD)

# The Matrix M is our Augmented Matrix for the Coefficients
M = []

# Here creates the two equations from the boundary conditions
for i in range(o):

	R=[]

	for j in range(N):

		R.append(pe(X[i],TDs[0][j]))


	M.append(R)

DX = (X[1]-X[0])/(N-o+1)
x = [DX*(i+1) for i in range(N-o)]

for i in range(N-o):
	R=[]
	for j in range(N):

		S = 0

		for k in range(len(a)):

			S += a[k]*pe(x[i],TDs[k][j])

		R.append(S)


	M.append(R)


for i in range(N-o):
	Y.append(f(x[i]))


M = np.array(M)
B = np.array(Y)

a = np.linalg.solve(M,B)

def af(x,a,C,N):
	S = 0
	for i in range(len(C)):
		S += a[i]*pe(x,C[i])

	return S



# Exact Plot 

def ef(x):
	e= math.e
	return (e**(2-x)-e**x)/(1-e**2)


points = 100
dx = (math.pi/2)/points
xlist = [dx*i for i in range(points+1)]
ylist = [ef(x) for x in xlist]
plt.plot(xlist,ylist, "-r")

# Our Approximation
yalist = [ af(x,a,TDs[0],N) for x in xlist]
plt.plot(xlist,yalist, "-g")
plt.show()

