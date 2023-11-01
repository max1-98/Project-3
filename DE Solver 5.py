# In this we create a callable function for evaluating linear ode's. Later we will modify this to add other basis functions, see DE Solver 6 when it's released.

from main import *
import mpmath
from mpmath import *
import math
import matplotlib.pyplot as plt

mpmath.mp.dps = 100
mpmath.mp.prec = 100

def DESolver(g,f,D,XY=[],DXDY=[], basis="c", N=20):

	# g = [g0,g1,g2,g3...]
	# f is what our ODE equals
	# D is our domain
	# XY is our initial coordinates
	# DXDY is our initial derivatives
	# basis represents what basis we will use, c means chebyshev

	# Order of the ODE
	n = len(g)-1

	# We start off by creating the Collocation points. 
	if basis == "c":
		# Collocation Points - Chebyshev
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


z = [mpf('-1.9446922743316067834825200168063e-62'), mpf('1.0'), mpf('-1.1739941671472299541508282885502e-30'), mpf('0.16691917482546069408301936969609'), mpf('-5.6584563617978052583906368581801e-29'), mpf('0.063972357794579552221197268497175'), mpf('1.0777456002327497351730488701764e-27'), mpf('0.2097373562230419938187855819954'), mpf('-7.9727737797990398340953597390266e-27'), mpf('-1.166908413765676939635859686587'), mpf('3.0235500782929361651708334096413e-26'), mpf('4.8369087620907615600665356528924'), mpf('-6.4423398376595195003165747179841e-26'), mpf('-11.306030908666381886967983165895'), mpf('7.8168056790060963502029830562517e-26'), mpf('15.553063360958604696794166007165'), mpf('-5.0528648615337027692563689553751e-26'), mpf('-11.562593179920282103494633395628'), mpf('1.3525685691207954964604042024982e-26'), mpf('3.6578645623387100944863035925434')]


def g0(x):
	return -mpmath.sin(polyeval(z,x))*polyeval(polydiff(z),x)

def g1(x):
	return mpmath.cos(polyeval(z,x))




# What our ODE equals
def f(x):

	return 1-mpmath.sin(polyeval(z,x))*polyeval(polydiff(z),x)*polyeval(z,x)
D = [-0.9,0.9]
g = [g0,g1]
solution = DESolver(g,f,D=D,XY=[[0],[0]])
points = 100
dx = (D[1]-D[0])/(points-1)

xlist = [D[0]+i*dx for i in range(points)]
def ef(x):
	return math.asin(x)

eylist = [ef(x) for x in xlist]
plt.plot(xlist,eylist,"-r")



aylist = [polyeval(solution,x) for x in xlist]
plt.plot(xlist,aylist,"-g")
plt.show()