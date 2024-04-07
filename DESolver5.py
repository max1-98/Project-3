# In this we create a callable function for evaluating linear ode's. Later we will modify this to add other basis functions, see DE Solver 6 when it's released.

from main import *
import mpmath
from mpmath import *
import math
import matplotlib.pyplot as plt
import latex
import numpy as np

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
	elif basis == "t":
		# Trigonometric
		points1 = ChebyshevNodes(m)

		# Transformed Points
		points2 = []

		for point in points1:
			points2.append((point*(b-a)/2+(b+a)/2))
	elif basis == "l":
		j = 1
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
	elif basis == "t":
		# Trigonometric basis the basis functions are very predictable so we will just form them as needed
		if m % 2 == 1:
			c = [0]*int((m-1)/2+1)
			s = [0]*int((m-1)/2)

		polys = []
		for i in range(m):
			c1 = c.copy()
			s1 = s.copy()
			if i < int((m+1)/2):
				c1[i] = 1
			else:
				s1[int(i-(m+1)/2)] = 1
			polys.append([c1,s1])
	elif basis == "l":
		return LaguerreGen(m)


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

	# We start off by creating the Collocation points and Basis Functions
	if basis == "c":
		# Chebyshev 
		points2 = points(N-n,a,b,"c")
		polys = polygen(N,a,b,"c")
		
	elif basis == "ec":
		# Even Chebyshev
		points2 = points(N-n,a,b,"ec")
		polys = polygen(N,a,b,"ec")
	elif basis == "t":
		points2 = points(N-n,a,b,"t")
		polys = polygen(N,a,b,"t")
	elif basis == "l":
		polys = polygen(N,a,b,"l")
		p = polys[N-2].copy()
		points2 = rootfind(p,NR=50)


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
	elif basis == "t":
		# Now we will repeat the above except for our Trignometric Polynomials. For simplicity in coding we will restrict the value
		# of N to being odd
		if N % 2 == 0:
			return "Error"

		# Initialise M and B
		M = []
		B = []
		if XY != []:
			for i in range(len(XY[0])):

				R=[]

				for j in range(int((N+1)/2)):
					R.append(cos(XY[0][i]*j))

				for j in range(int((N-1)/2)):
					R.append(sin(XY[0][i]*(j+1)))

				B.append(XY[1][i])
				M.append(R)


		# Derivatives
		if DXDY != []:
			for i in range(len(DXDY[0])):

				R = []

				for j in range(int((N+1)/2)):
					R.append(-i*sin(DXDY[0][i]))

				for j in range(int((N-1)/2)):
					R.append((i+1)*cos((i+1)*DXDY[0][i]))


				B.append(DXDY[1][i])
				M.append(R)

		for point in points2:
			R = []

			for poly in polys:
				s = 0
				for i in range(n+1):
					s += g[i](point)*trigeval(trigdiffn(poly.copy(),i),point)
				R.append(s)

			B.append(f(point))
			M.append(R)



	M = mpmath.matrix(M)
	B = mpmath.matrix(B)


	A = mpmath.lu_solve(M,B)

	if basis != "t":
		solution = []
		for i in range(N):
			solution = polyadd(solution,polysmult(polys[i],A[i]))
	else:
		solution = [A[0:int((N+1)/2)],A[int((N+1)/2):N]]

	return solution

