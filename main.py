# This document will contain all of my functions that I will call upon later

import math
import numpy
import mpmath
mpmath.mp.dps = 100
mpmath.mp.prec = 100

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
		CN.append(mpmath.cos(((1+2*i)*mpmath.pi)/(2*N)))

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


# Here we create our Non-Linear Representations along with their Operations so we can generate our system of equations
# from an ODE. We'll also introduce a Jacobian Function
# See my Notes for how each of these operators is meant to work mathematically

def dim(x):

	z = x.copy()
	n = 0
	while not isinstance(z, int) and not isinstance(z,int):
		z = z[-1]
		n += 1


	return n

def add(x,y):
	n = len(x)
	m = len(y)

	if n > m:
		m = n

	for i in range(m):

		# Then we just term wise add

		for j in range(len(x[i])):

			if len(x[i]) != len(y[i]):
				return "Size-Error"
			if i < 2:
				x[i][j] = x[i][j] + y[i][j]
			elif i < 3:
				for k in range(len(x[i])):
					x[i][j][k] = x[i][j][k] + y[i][j][k]

	return x

def smult(x,s):

	n = len(x)

	for i in range(n):

		# Then we just term wise add

		for j in range(len(x[i])):

			if i < 2:
				x[i][j] = s*x[i][j]
			elif i < 3:
				for k in range(len(x[i])):
					x[i][j][k] = s*x[i][j][k]
	return x


def mult(x,s):
	n = dim(x)

	if n == 1:
		for i in range(len(x)):
			x[i] = x[i]*s
	elif n == 2:
		for i in range(len(x)):
			for j in range(len(x)):
				x[i][j] = x[i][j]*s

	return x

def mcross(x,y):
	# This is that product rule described in my Notes except just for when they are singular ones.
	# i and j are the sizes of x and y respectively

	# First we need to determine the size of x and y
	n = dim(x)
	k = dim(y)

	for i in range(n+k):
		m = [0]*len(y)
		for j in range(i):
			M = []
			for k in range(len(y)):
				M.append(m.copy())

			m = M.copy()
	if k+n == 0:
		return [x[0]*y[0]]

	if k+n == 1:
		if n > k:
			print(x)
			print(y)
			return mult(x,y[0])
		else:
			return mult(y,x[0])

	if k+n == 2:
		for i in range(len(x)):
			for j in range(len(x)):
				m[i][j] = x[i]*y[j]

	return m


def cross(x,y):
	# Our hardest function. This will deal with multiplying full representations of what we saw out

	# Our new size
	n = len(x)+len(y)-2


	# So here we initialize our result. We create it as just zeros. 

	z = [[0]]
	for i in range(n):
		m = [0]*n
		for j in range(i):
			M = []
			for k in range(n):
				M.append(m.copy())

			m = M.copy()
		z.append(m.copy())

	for xs in x:
		for ys in y:
			xs2 = xs.copy()
			ys2 = ys.copy()
			z = add(mcross(xs2,ys2),z)

	return z 

#print(cross([[0],[0,1],[[1,1],[1,1]]],[[0],[0,1],[[1,1],[1,1]]]))


def integrate(f,a,b,N=1000):

	dx = (b-a)/N
	S = f(a)+f(b)

	for i in range(1,N):
		S += 2*f(a+i*dx)


	return 1/2*dx*S
