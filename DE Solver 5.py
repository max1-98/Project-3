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


z =[mpf('1.0000000000000000000000000000032'), mpf('-2.6116399978094670744644415892405e-31'), mpf('-0.16666518468218911239248943440958'), mpf('-0.000047835981503579706215130099133959'), mpf('0.016969667291367152098645947876485'), mpf('-0.00082768565738339960510628876686211'), mpf('-0.00021315033106398116541657633767751'), mpf('-0.0011422111214554475206250571831621'), mpf('0.00079340538924191988255744007993649'), mpf('-0.00025719796029884002148253632830167'), mpf('0.000049292913255131130572971688965504'), mpf('-0.0000055123584567614088154532742009583'), mpf('0.00000021699557462006678757781513248084'), mpf('0.000000035770139645078675023748111660065'), mpf('-0.0000000071864691730238485756714847841217'), mpf('6.3167756019371508807520735137327e-10'), mpf('-3.1661312697384011647279135234403e-11'), mpf('8.5362520164817496448993247615638e-13'), mpf('-8.0781727545594862508031767606753e-15'), mpf('-6.7176119461117863016612977665737e-17')]
plt.title(r"$y_4$")
def g0(x):
	return x*2*polyeval(z,x)**1

def g1(x):
	return 2

def g2(x):
	return x

# What our ODE equals
def f(x):

	return x*1*polyeval(z,x)**2


D = [0,10]
g = [g0,g1,g2]
solution = DESolver(g,f,D=D,XY=[[0],[1]],DXDY=[[0],[0]],basis="c",N=20)
print(solution)

points = 100
dx = (D[1]-D[0])/(points-1)

xlist = [D[0]+i*dx for i in range(points)]



def ef(x):
	return polyeval([mpf('0.99999999999999999999999999999606'), mpf('-3.5883422671845345817255564032015e-31'), mpf('-0.16666518468218911239984003641286'), mpf('-0.000047835981503579468270797760193845'), mpf('0.016969667291367150569200473249448'), mpf('-0.00082768565738339530328224136872361'), mpf('-0.00021315033106398795342708813115583'), mpf('-0.0011422111214554407301651494157264'), mpf('0.0007934053892419152579808237879'), mpf('-0.0002571979602988377787236599104652'), mpf('0.000049292913255130333236865966524957'), mpf('-0.0000055123584567611972219209250143059'), mpf('0.00000021699557462002447178529679580262'), mpf('0.000000035770139645085064126338860815762'), mpf('-0.0000000071864691730245723630540727197657'), mpf('6.3167756019377562047074031157276e-10'), mpf('-3.1661312697387636088660845843879e-11'), mpf('8.5362520164832189276688000528998e-13'), mpf('-8.0781727545630977022428348259323e-15'), mpf('-6.7176119461077230837351345430033e-17')],x)

eylist = [ef(x) for x in xlist]
plt.plot(xlist,eylist,"-r",label=r"$y*$")



aylist = [polyeval(solution,x) for x in xlist]
plt.plot(xlist,aylist,"--g",label=r"$y=y_4$")

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
