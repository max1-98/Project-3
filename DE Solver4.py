

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
N = 50

# Coefficients of the terms in the ODE "Coeffunctions"
# From y, dy/dx, d^2y/dx^2...
z = [mpf('0.99999999999999999999999999999921'), mpf('1.9647051309752676916335479445083e-49'), mpf('0.49999999999999999884135053710615'), mpf('1.3762049078006378300066801029233e-30'), mpf('-0.12499999999999976743322573415228'), mpf('-3.6548557930630364442700693102597e-28'), mpf('0.062499999999981338178776596500379'), mpf('3.2376868499803112024689413518296e-26'), mpf('-0.039062499999200509200039312784485'), mpf('-1.4620561223847287700247488583099e-24'), mpf('0.027343749978828934752406203223289'), mpf('3.9423585806660005344278036806864e-23'), mpf('-0.020507812121494563748732470968982'), mpf('-6.9322020543106262791838750687996e-22'), mpf('0.016113276405613427626666052622584'), mpf('8.4606532492836200626714587458595e-21'), mpf('-0.013091994756822020148274762240591'), mpf('-7.5031783968656868302225460799461e-20'), mpf('0.010909694434958476982260001305353'), mpf('5.0014457604128844126020601918176e-19'), mpf('-0.00927156505376515934085055592754'), mpf('-2.5695580349723333498574389342565e-18'), mpf('0.007999858752562457069663937347276'), mpf('1.0365611025601279310554103476002e-17'), mpf('-0.0069735778000362199054845917402771'), mpf('-3.3278031501862458204801577778126e-17'), mpf('0.006093150654547874715770691061813'), mpf('8.5822201198476463121725456841214e-17'), mpf('-0.0052621862851085421691365811085257'), mpf('-1.7882146429818538079656387628628e-16'), mpf('0.0043934646005613559131250997521427'), mpf('3.0177097041283868881409961246017e-16'), mpf('-0.003442709812151591099671233212248'), mpf('-4.1208871224253784667308302538556e-16'), mpf('0.0024470404280883868179532156123677'), mpf('4.5331470015622872158033492958702e-16'), mpf('-0.0015223645278886370461926117514901'), mpf('-3.9815848298507139867409348647292e-16'), mpf('0.00079908953487350992454840145486511'), mpf('2.7518786677774123853564433007309e-16'), mpf('-0.00034004349339141970010946067722129'), mpf('-1.4628711322438159036119977309343e-16'), mpf('0.00011172211413379445018296139222167'), mpf('5.7692993304996739045945754448196e-17'), mpf('-0.000026434307142543192838598375051487'), mpf('-1.5892814027608192371664022605086e-17'), mpf('0.000003991678752886522507925728271197'), mpf('2.728671450809787585448783752665e-18'), mpf('-0.00000028805280642000963895130054472721'), mpf('-2.1975623073275430285429251060722e-19')]
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