# In this we create a callable function for evaluating linear ode's. Later we will modify this to add other basis functions, see DE Solver 6 when it's released.

from main import *
import mpmath
from mpmath import *
import math
import matplotlib.pyplot as plt
import latex

plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.grid()

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
		polys = polygen(N,a,b,"ec")





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


z = [mpf('1.0'), mpf('1.0746689089680776060052356919624e-31'), mpf('-0.16672020794374009368173069114232'), mpf('0.0017347969670028699952536258448428'), mpf('0.030476776882736255008442090940891'), mpf('0.031367983936293708132863028796137'), mpf('-0.059705443130453865666436250207982'), mpf('0.044138992383643465702358398679791'), mpf('-0.020878146958221387075262241860894'), mpf('0.0071705769629258347163718694948124'), mpf('-0.0018853707284993811178918233806489'), mpf('0.00038852921454656700668525793713618'), mpf('-0.000063189238100918680871747046644829'), mpf('0.0000080744578331789335381866759313996'), mpf('-0.00000079953271712729191308361756686267'), mpf('0.000000059920184102454497172442403736275'), mpf('-0.0000000032749187293368076399349849567233'), mpf('1.2279767047232116846275923165332e-10'), mpf('-2.8177206999453074313004284337679e-12'), mpf('2.978001446074926878769596543453e-14')]


plt.title(r"$y_6$")
def g0(x):
	return x*5*polyeval(z,x)**4

def g1(x):
	return 2

def g2(x):
	return x

# What our ODE equals
def f(x):

	return x*4*polyeval(z,x)**5


D = [0,10]
g = [g0,g1,g2]
solution = DESolver(g,f,D=D,XY=[[0],[1]],DXDY=[[0],[0]],basis="c",N=20)
print(solution)

points = 100
dx = (D[1]-D[0])/(points-1)

xlist = [D[0]+i*dx for i in range(points)]


def ef(x):
	return 1/(math.sqrt(1+x**2/3))

eylist = [ef(x) for x in xlist]
plt.plot(xlist,eylist,"-r",label=r"$y=\frac{1}{\sqrt{1+\frac{x^2}{3}}}$")



aylist = [polyeval(solution,x) for x in xlist]
plt.plot(xlist,aylist,"--g",label=r"$y=y_6$")

"""
polys=[[1],[-1,0,2],[1,0,-8,0,8]]
a = [    1.2173411583055409963875616632, 0.20842027928171569810339479336, -0.0089208790238252982841668698426]
solution1 = []

for i in range(3):
	solution1 = polyadd(solution1,polysmult(polys[i],a[i]))
print(solution1)
def f(x):

	return polyeval(solution1,x)
nrlist = []
for x in xlist:
	nrlist.append(f(x))
plt.plot(xlist,nrlist,"--b")
"""
plt.legend(loc="upper right")
plt.show()


# We can call the above function and it will plot our ODE 
