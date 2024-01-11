# This document will contain all of my functions that I will call upon later

import math
import numpy as np
import mpmath
from mpmath import *
import cmath
mpmath.mp.dps = 100
mpmath.mp.prec = 100

### Polynomial Functions ###

def bexp(a,n):
	# Expands (a+x)^n and leaves it as a polynomial

	p = []

	for i in range(n+1):
		p.append(a**(n-i)*math.factorial(n)/(math.factorial(n-i)*math.factorial(i)))

	return p

# This function adds two polynomials
def polyadd(p1,p2):
	l1 = len(p1)
	l2 = len(p2)

	if l1 > l2:
		for i in range(l2):
			p1[i] = p1[i]+p2[i]
	else:
		for i in range(l1):
			p2[i] = p1[i]+p2[i]	

		p1 = p2.copy()

	return p1


# This function differentiates a polynomials
def polydiff(F):

	# F is a list of coefficient of the polynomial a0+a1x^1...
	l = len(F)
	f = [0]*len(F)

	for i in range(l-1):

		f[i] = (i+1)*F[i+1]

	return f

def polydiffn(F,n):

	for i in range(n):
		F = polydiff(F)

	return F

# Polynomial evaluator 
def polyeval(a,x):
	# Polynomial Eval
	n = len(a)

	S = 0 
	for i in range(n):

		S += a[i]*x**i

	return S

# Integral function, f is the list of coefficients of our polynomial
def polyint(f):

	F = (len(f)+1)*[0]

	for i in range(len(f)):
		F[i+1] = f[i]/(i+1)

	return F


# This function multiplies two polynomials together
def polymult(p1,p2):

	l1 = len(p1)
	l2 = len(p2)
	p3 = (len(p1)+len(p2)-1)*[0]

	for i in range(l1):
		for j in range(l2):

			p3[i+j] += p1[i]*p2[j]

	return p3

# This will multiply polynomials together and return the result 
def polysmult(p1, a):

	for i in range(len(p1)):
		p1[i] = a*p1[i]

	return p1

def polyip(p1,p2,W="L"):

	if W == "L":

		p0 = polyint(polymult(p1,p2))
		
		return polyeval(p0,1)-polyeval(p0,-1)

# Transformations
def polyxstret(p,a):
	# Computes the calculation p(x) -> p(x/a)
	o = len(p)

	for i in range(o):
		p[i] = p[i]/a**i

	return p

def polyystret(p,a):
	# Computes the calculation p(x) -> ap(x)
	o = len(p)

	for i in range(o):
		p[i] = a*p[i]

	return p

def polyxtrans(p,a):
	# p(x) -> p(x-a)
	o = len(p)

	p1 = []

	for i in range(o):
		p2 = bexp(-a,i)
		p2 = polysmult(p2,p[i])
		p1 = polyadd(p1,p2)
	
	return p1

def polyytrans(p,a):
	# p(x) -> p(x) + a

	p[0] += a

	return p

### Orthogonal Polymials ###


# Chebyshev Polynomials
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

# Chebyshev Nodes 
def ChebyshevNodes(N):

	CN = []
	for i in range(N):
		CN.append(mpmath.cos(((1+2*i)*mpmath.pi)/(2*N)))

	return CN

# Laguerre Polynomial
def LaguerreGen(N):

	L = [[1],[1,-1]]

	for i in range(N-2):
		L.append(polysmult(polyadd(polymult([2*i+3,-1],L[i+1].copy()),polysmult(L[i].copy(),-(i+1))),1/(i+2)))

	return L

### Root Finding 

# High Accuracy Root Finding For a Polynomial
def rootfind(p,NR=40):
	# Our guesses will come from using the np.roots()
	# We'll use Numpy as we aren't too bothered about the accuracy
	# p = [p0,p1...pn]
	p1 = p.copy()

	p1.reverse()
	poly = np.poly1d(p1)
	roots = poly.roots

	roots1 = []
	for root in roots:
		roots1.append(root.real)

	roots1.sort()

	# Next for additional accuracy we can use the Newton-Raphson method on each root
	dp = polydiff(p.copy())
	for i in range(len(p)-1):

		x = mpf(roots1[i])
		for j in range(NR):
			x = x-polyeval(p,x)/polyeval(dp,x)

		roots1[i] = x

	return roots1


#--- The Following Function is Unstable Numerically ---#
def rootsqr(p,QR=100):
	# First we convert p to a monic polynomial
	p = polysmult(p,mpf(1/p[-1]))
	N = QR
	DIM = len(p)-1

	# Then we form our Companion Matrix
	ROW = [0]*DIM

	A = []
	for i in range(DIM):
		A.append(ROW.copy())

	for j in range(DIM-1):
		A[j+1][j] = 1

	A = np.array(A,dtype=object)

	for i in range(DIM):
		A[i,DIM-1] = -p[i]

	# Now we use the QR Algorithm to find an approximation of the roots
	for k in range(N):
		# Next we need to find out Orthogonal Matrix. 
		Q = A.copy()

		qm = []

		for i in range(DIM):

			# ith row of A
			a = A[:,i].copy()

			# ith row of A, used to initialize ith row of q
			qt = A[:,i].copy()

			# Applying G-S Algorithm
			for j in range(i):
				qt += -np.dot(a,Q[:,j].copy())*Q[:,j].copy()

			# Normalising qt
			qm.append(np.dot(qt,qt))
			q = (1/(qm[i])**(1/2))*qt

			# Setting the ith row of Q to be q
			Q[:,i] = q

		# Now working on R
		R = np.zeros([DIM,DIM])

		for i in range(DIM):
			for j in range(DIM):
				if j >= i:	
					if i == j:
						R[i,j] = (qm[i])**(1/2)
					else:
						R[i,j] = np.dot(A[:,j],Q[:,i])

		# Now we find the next term by computing RQ
		A = np.matmul(R,Q).copy()

	# The approximate roots are the diagonal elements of A	
	roots = []
	for i in range(DIM):
		roots.append(A[i,i])
#-------------------------------------------------------#


### Numerical Integration ###

# Trapezium Rule
def integrate(f,a,b,N=1000):

	dx = (b-a)/N
	S = f(a)+f(b)

	for i in range(1,N):
		S += 2*f(a+i*dx)


	return 1/2*dx*S

### Numerical Differentiation ### 
def CSD(r,x,a=0,w=1,d=1,h=10**(-8)):


	# f is the function, x is the point we are calculating the derivative at and d specifies which derivative
	def f(x,r,a=0,w=1):

		if x == a+r*w:
			return 1
		else:
			return sin(pi/w*(x-a-r*w))/(pi/w *(x-a-r*w))

	if d == 1:
		z = complex(x,h)
		return f(x,r,a,w).imag/h
	elif d == 2:
		z = complex(x,h)
		return -2*(f(z,r,a,w).real-f(x,r,a,w))/h**2
	elif d == 3:
		z = complex(x,h)
		y = complex(x,2*h)
		return 1/h**3*(2*f(z,r,a,w).imag-f(y,r,a,w).imag)

	return "Please select an integer d between 0<d<4"


### Trigonometric Functions ###

def trigeval(p,x):
	# In this case our list is a nested list: [[an],[bn]]
	# With an the coefficients of the cosines and bn the coefficients of the sines
	# The convention is b0=0 so we will ignore this to avoid invertibility problems later
	# [[a0,a1,a2,a3...],[b1,b2,b3...]]

	s = 0

	for i in range(len(p[0])):
		s += p[0][i]*cos(i*x)

	for i in range(len(p[1])):
		s += p[1][i]*sin((i+1)*x)

	return s

def trigdiff(p):
	# In this case our lists get swapped in order and multiplied by the i value
	pc = []

	# Including the cos(0x) term which goes to 0
	pc0 = [0]

	# Using Derivative of sin(nx) -> ncos(nx)
	for i in range(len(p[1])):
		pc0.append(p[1][i]*(i+1))

	pc1 = []
	# Using Derivative of cos(nx) -> -nsin(nx)
	for i in range(1,len(p[0])):
		pc1.append(-p[0][i]*i)

	# The Sines become the cosines and vice versa
	pc.append(pc0)
	pc.append(pc1)

	return pc

def trigdiffn(p,n):

	pc = p
	# Note the for loop won't do anything if n=0, corresponding to the 0th derivative ie. just our function
	for i in range(n):
		pc = trigdiff(pc)

	return pc

def trigadd(p1,p2):
	c1 = p1[0]
	s1 = p1[1]
	c2 = p2[0]
	s2 = p2[1]

	if len(c1) > len(c2):
		for i in range(len(c1)):
			c3 = c1.copy()
			c3[i] = c1[i] + c2[i]
	else:
		for i in range(len(c2)):
			c3 = c2.copy()
			c3[i] = c1[i] + c2[i]


	if len(s1) > len(s2):
		s3 = s1.copy()
		for i in range(len(s1)):
			s3[i] = s1[i] + s2[i]
	else: 
		s3 = s2.copy()
		for i in range(len(s2)):
			s3[i] = s1[i] + s2[i]
	return [c3,s3]

def trigmult(p,m):
	p1 = p.copy()
	for i in range(2):
		for j in range(len(p[i])):
			p1[i][j] = m*p[i][j]

	return p1

