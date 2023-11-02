# In this we create a callable function for evaluating linear ode's. Later we will modify this to add other basis functions, see DE Solver 6 when it's released.

from main import *
import mpmath
from mpmath import *
import math
import matplotlib.pyplot as plt

mpmath.mp.dps = 100
mpmath.mp.prec = 100

def points(m,a,b,basis="c"):
	# Collocation Points - Chebyshev
	if basis == "c":
		# Chebyshev

		points1 = ChebyshevNodes(m)

		# Transformed Points
		points2 = []

		for point in points1:
			points2.append((point*(b-a)/2+(b+a)/2))
	elif basis == "ec":
		# Even Chebyshev

		points1 = ChebyshevNodes(2*m)

		# Transformed Points
		points2 = []

		for point in points1:
			if point > 0:
				points2.append((point*(b-a)/2+(b+a)/2))

	return points2

def polygen(m,a,b,basis="c"):

	if basis =="c":
		# Chebyshev

		polysNT = ChebyshevGen(m)

			#Creating Our Translated Versions
		polys = []
		for i in range(m):
			p1 = polysNT[i]
			p1 = polyxstret(p1,(b-a)/2)
			p1 = polyxtrans(p1,(b+a)/2)
				
			polys.append(p1)
	elif basis == "ec":
		# Even Chebyshev
		polysNTOE = ChebyshevGen(2*m)

		# Now we remove all the odd ones
		i = 0
		polysNT = []
		for poly in polysNTOE:

			if i % 2 == 0:
				polysNT.append(poly)

			i += 1


		#Creating Our Translated Versions
		polys = []
		for i in range(m):
			p1 = polysNT[i]
			p1 = polyxstret(p1,(b-a)/2)
			p1 = polyxtrans(p1,(b+a)/2)
				
			polys.append(p1)

	return polys


def DESolver(g,f,D,XY=[],DXDY=[], basis="c", N=20):

	# g = [g0,g1,g2,g3...]
	# f is what our ODE equals
	# D is our domain
	# XY is our initial coordinates
	# DXDY is our initial derivatives
	# basis represents what basis we will use, c means chebyshev

	# Order of the ODE
	n = len(g)-1
	a = D[0]
	b = D[1]

	# We start off by creating the Collocation points. 
	if basis == "c":
		# Chebyshev 
		points2 = points(N-n,a,b,"c")
		polys = polygen(N,a,b,"c")
		
	elif basis == "ec":
		# Even Chebyshev
		points2 = points(N-n,a,b,"ec")
		print(points2)
		polys = polygen(N,a,b,"ec")
		print(polys)



	# So polys are our transformed Chebyshev Polynomials

	### Next let's create the N equations

	# This works so long as it's not a trigonometric polynomial basis
	if basis != "t":
		# Mx=B where x is the vector representing our coefficients a0,a1,a2...
		M = []
		B = []

		# We'll start with the Initial Conditions

		# Coordinates
		if XY != []:
			for i in range(len(XY[0])):

				R=[]

				for j in range(N):

					R.append(polyeval(polys[j],XY[0][i]))
				B.append(XY[1][i])
				M.append(R)


		# Derivatives
		if DXDY != []:
			for i in range(len(DXDY[0])):

				R = []
				polysd1 = polys.copy()
				polysd = []

				for poly in polysd1:
					polysd.append(polydiff(poly))

				for j in range(N):
					R.append(polyeval(polysd[j],DXDY[0][i]))


				B.append(DXDY[1][i])

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
					

					S += g[k](x)*polyeval(p2,x)

				R.append(S)

			M.append(R)
			B.append(f(x))


	M = mpmath.matrix(M)
	B = mpmath.matrix(B)
	A = mpmath.lu_solve(M,B)


	solution = []
	for i in range(N):
		solution = polyadd(solution,polysmult(polys[i],A[i]))

	return solution


z = [mpf('1.0'), mpf('0.0'), mpf('0.49999999999999999900863376932315'), mpf('0.0'), mpf('-0.12499999999999964162110761058566'), mpf('0.0'), mpf('0.062499999999971190153811881013597'), mpf('0.0'), mpf('-0.031249999998963480631219292746278'), mpf('0.0'), mpf('0.015624999978922638322370276372909'), mpf('0.0'), mpf('-0.0078124997267616825809087939781773'), mpf('0.0'), mpf('0.0039062475707294674989700189517534'), mpf('0.0'), mpf('-0.0019531094562936475527404081076312'), mpf('0.0'), mpf('0.00097648849599602535579593078509154'), mpf('0.0'), mpf('-0.00048801282507110525479452705899763'), mpf('0.0'), mpf('0.00024338633196584852509425904915495'), mpf('0.0'), mpf('-0.00012040871337425147468416429225482'), mpf('0.0'), mpf('0.000058142359188253027723163233567748'), mpf('0.0'), mpf('-0.000026514302695524346876603697214962'), mpf('0.0'), mpf('0.000010832853251309743207544779043098'), mpf('0.0'), mpf('-0.0000036922491001606837038059953517545'), mpf('0.0'), mpf('0.0000009562264724059215301192160864378'), mpf('0.0'), mpf('-0.00000016350280470190300817067341614658'), mpf('0.0'), mpf('0.000000013625233725158573081076786296374'), mpf('0.0')]
def g0(x):
	return polyeval(polydiff(z),mpf(x))

def g1(x):
	return polyeval(z,mpf(x))


# What our ODE equals
def f(x):

	return mpf(x)+polyeval(z,mpf(x))*polyeval(polydiff(z),mpf(x))


D = [-1,1]
g = [g0,g1]
solution = DESolver(g,f,D=D,XY=[[0],[1]],basis="ec",N=3)

points = 100
dx = (D[1]-D[0])/(points-1)

xlist = [D[0]+i*dx for i in range(points)]
def ef(x):
	return math.sqrt(1+x**2)

eylist = [ef(x) for x in xlist]
plt.plot(xlist,eylist,"-r")



aylist = [polyeval(solution,x) for x in xlist]
plt.plot(xlist,aylist,"-g")
plt.show()


# We can call the above function and it will plot our ODE 