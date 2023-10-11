# This document will contain all of my functions that I will call upon later

import math
import numpy

###*** Below are all of the polynomial functions we will need ***###

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

# Orthogonal Polymials 

# Chebyshev Polynomials Generator
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

# This finds the Nth set of Chebyshev Nodes 
def ChebyshevNodes(N):

	CN = []
	for i in range(N):
		CN.append(math.cos((math.pi+2*math.pi*i)/(2*N)))

	return CN


# This creates the ith Lagrange Interpolating polynomial of a set of x values called x
def lagpoly(x, i):

	N = len(x)

	# Here we create our factor. In this case we have calculate P but later we will multiply by 1/P
	P = 1
	for j in range(N):
		if j != i:
			P *= (x[i]-x[j])

	# Initialise our polynomial as p(x) = 1
	p1 = [1]

	# Now we multiply by each factor
	for j in range(N):
		if j != i:
			p1 = polymult(p1,[-x[j],1])

	return polysmult(p1, 1/P)

def legpoly(N):

	# This will generate the first N normalised legendre polynomials N > 2
	P = [[1/math.sqrt(2)]]

	# We'll be doing this by Gram Schmidt
	for n in range(1,N):

		# Stage 1 
		p1 = [0]*(n)
		p1.append(1)


		# Stage 2 Creating our orthogonal polynomial
		p3 = p1
		for j in range(n):

			p2 = P[j].copy()
			p2 = polysmult(p2,-polyip(p2,p3))
			p1 = polyadd(p1,p2)


		# Stage 3 Normalising our polynomial
		P.append( polysmult(p1,1/math.sqrt(polyip(p1,p1))) )

	return P


def gaussint(f,a,b, N=10, nodes="C"):

	# This function will approximate the integral between a and b of f(x) using Gaussian Integration


	### N is the number of nodes

	if nodes == "C":
		# Chebyshev

		# NOTE Chebyshev GI Weights are always pi/n
		CN = ChebyshevNodes(N)

		# We need to translate the CN first. 
		S = 0

		for i in range(N):

			S += f((b-a)/2*CN[i]+(b+a)/2)*math.sqrt(1-CN[i]**2)

		return (b-a)/2*math.pi/N*S
	elif nodes == "L":
		# Lobatto
		LN = legpoly(N)
	else:
		return 0


