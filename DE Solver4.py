

# In this again we will solve an N dimensional linear DE. 


from main import *
import math
import numpy as np
import matplotlib.pyplot as plt
import mpmath
from mpmath import *

mpmath.mp.dps = 100
mpmath.mp.prec = 100


### Inputs from the User

# Accuracy (number of Orthogonal Polys)
N = 20

# Coefficients of the terms in the ODE "Coeffunctions"
# From y, dy/dx, d^2y/dx^2...
z = [mpf('1.0'), mpf('-7.5964541966078389979785938156496e-65'), mpf('0.5'), mpf('6.343839899453709649337223491084e-34'), mpf('-1.3184491312894828952604967780596e-32'), mpf('-2.2262841024800207473019877057249e-32'), mpf('1.4441277150173815430340824006955e-31'), mpf('2.6311439424014895920469823713884e-31'), mpf('-7.7377855768554919190097963986379e-31'), mpf('-1.4787079897576693696382051402784e-30'), mpf('2.2796588964003716362906393884204e-30'), mpf('4.5240751362014822131501800383165e-30'), mpf('-3.8783095199612239091355305752283e-30'), mpf('-7.957031190709998956810081033472e-30'), mpf('3.8015281196433160340327112903093e-30'), mpf('8.0363161052596793013838729633121e-30'), mpf('-1.9955092714426741638710618253969e-30'), mpf('-4.3316969961217474361602260341421e-30'), mpf('4.3475942422129448953908708629608e-31'), mpf('9.6563660226763096263384428635804e-31')]



def g0(x):
	return polyeval(polydiff(z),x)

def g1(x):
	return polyeval(z,x)
def g2(x):
	return (1-x**2)

def g3(x):

	return x**4


# What our ODE equals
def f(x):

	return polyeval(z,x)*polyeval(polydiff(z),x)+x

# Initial Conditions

X = [0]
Y = [1]

X1 = []
DY = []

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

# Coordinates
for i in range(len(X)):

	R=[]

	for j in range(N):

		R.append(polyeval(polys[j],X[i]))
	B.append(Y[i])
	M.append(R)


# Derivatives
for i in range(len(X1)):

	R = []
	polysd1 = polys.copy()
	polysd = []

	for poly in polysd1:
		polysd.append(polydiff(poly))

	for j in range(N):
		R.append(polyeval(polysd[j],X1[i]))


	B.append(DY[i])

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


# Change this to solve using mpf 
M = mpmath.matrix(M)
B = mpmath.matrix(B)
A = mpmath.lu_solve(M,B)


solution = []
for i in range(N):
	solution = polyadd(solution,polysmult(polys[i],A[i]))

points = 100
dx = (b-a)/(points-1)

xlist = [a+i*dx for i in range(points)]

"""
# Solution to dy/dx=y y(0)=1
def ef(x):

	return math.e**(-x)

eylist = [ef(x) for x in xlist]

"""


def ef(x):
	return math.sqrt(1+x**2)

eylist = [ef(x) for x in xlist]
plt.plot(xlist,eylist,"-r")

def nf(x):
	return 1+0.4797*x**2-0.06334*x**4

nylist = [nf(x) for x in xlist]
plt.plot(xlist,nylist,"-b")

aylist = [polyeval(solution,x) for x in xlist]

"""
def j(x,M=100):

	

	return -2*x**2+1

eylist = [j(x) for x in xlist]

plt.plot(xlist,eylist,"-r")
"""
plt.plot(xlist,aylist,"-g")
plt.show()
print(solution)
# p1 = [0.0, 1.850815717680926, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# p2 = [3.944304526105059e-31, 1.0000000000000009, -6.918772340517001e-14, 24.389568917486372, 4.2771433271321947e-13, -167.31904436954835, -7.22007280590506e-13, 294.23019433219423, 3.478214081967499e-13, -144.09649689397693]
# p3 = [0.0, 1.0, -2.4067833853528176e-13, -77.22848071406378, 1.353899045419129e-12, 539.065483032338, -2.3203773526874835e-12, -1063.0477066692667, 1.178724417193309e-12, 627.2597267024689]