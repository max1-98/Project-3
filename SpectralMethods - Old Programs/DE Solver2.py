
# in this we are going to solve y''+y=x^2 with y(0)=-1 and y(pi/2)=(pi/2)^2-2
# but in this example we will do so using the Chebyshev polynomials and the collocation points will be the Chebyshev nodes

import numpy as np
import math
import matplotlib.pyplot as plt

# Exact Plot 

def ef(x):
	return math.cos(x)+x**2-2


points = 100
dx = (math.pi/2)/points
xlist = [dx*i for i in range(points+1)]
ylist = [ef(x) for x in xlist]
plt.plot(xlist,ylist, "-r")



# Functions

def deriv(x):

	# x is a list of coefficient of the polynomial a0+a1x^1...
	l = len(x)
	y = [0]*len(x)

	for i in range(l-1):

		y[i] = (i+1)*x[i+1]

	return y

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

def pe(x,a):
	# Polynomial Eval
	n = len(a)

	S = 0 
	for i in range(n):

		S += a[i]*x**i

	return S


def SDChebyshev(C):
	# Create the second derivatives of the Chebyshev Polynomials
	C1 = []
	for T in C:
		C1.append(deriv(deriv(T)))

	return C1

# Our Chebyshev Spectral Approximation 

N = 10
C = ChebyshevGen(N)
C2 = SDChebyshev(C)

M = []
# X values of the initial Conditions, This forms our initial conditions
X = [0, math.pi/2]
Y = [-1, (math.pi/2)**2-2]
for i in range(2):
	R=[]
	for j in range(N):
		R.append(pe(X[i],C[j]))
	M.append(R)

DX = (X[1]-X[0])/(N-1)
x = [DX*(i+1) for i in range(N-2)]


for i in range(N-2):
	R=[]
	for j in range(N):
		R.append(pe(x[i],C2[j])+pe(x[i],C[j]))


	M.append(R)


for i in range(N-2):
	Y.append(x[i]**2)


M = np.array(M)
B = np.array(Y)

a = np.linalg.solve(M,B)

def af(x,a,C,N):
	S = 0
	for i in range(len(C)):
		S += a[i]*pe(x,C[i])

	return S


yalist = [ af(x,a,C,N) for x in xlist]
plt.plot(xlist,yalist, "-g")
plt.show()
