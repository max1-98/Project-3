# In this we will solve a simple set of equations, see my Notes for the full equations. We will copy and expand upon
# 3DegNewtonMethod2 
# Result:

# Error = 0.0000245848
# Interpolating Polynomial = [mpf('0.99999999999999999999999999999921'), mpf('0.0'), mpf('0.49973666857582575077137573566055'), mpf('0.0'), mpf('-0.12062258563421877376961495106674'), mpf('0.0'), mpf('0.045535878902564385505941386762077'), mpf('0.0'), mpf('-0.010425003185446731647617695197307'), mpf('0.0')]

from mpmath import *
from main import * 
mpmath.mp.dps = 100
mpmath.mp.prec = 100
import matplotlib.pyplot as plt

# We'll set this as our number of iterations
N = 50


# Need the 49 Cheybshev nodes. Ie. the positive roots of the 50th Chebyshev polynomial
# Need 50 Chebyshev polynomials 
x = [mpf('0.98078528040323044912618223613438'), mpf('0.83146961230254523707878837761816'), mpf('0.55557023301960222474283081394853'), mpf('0.19509032201612826784828486847702')]

b = [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [-1.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, -8.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0], [-1.0, 0.0, 18.0, 0.0, -48.0, 0.0, 32.0, 0.0, 0.0, 0.0], [1.0, 0.0, -32.0, 0.0, 160.0, 0.0, -256.0, 0.0, 128.0, 0.0]]

# Creates our f
def f(a,x,b):
	# Each row of f
	S1 = 0
	for i in range(5):
		S1 += a[i]*polyeval(b[i],0)

	out = [S1-1]

	for i in range(len(x)):

		# Initialises the row as 0
		S = 0

		for j in range(len(a)):
			#Sums the terms where the index is the same a0a0 
			S += (a[j]**2)*polyeval(polydiff(b[j]),x[i])*polyeval(b[j],x[i])

			# Sums the terms where they are different
			for k in range(len(a)):
				if j != k:
					S += a[j]*a[k]*(polyeval(polydiff(b[j]),x[i])*polyeval(b[k],x[i]))

		out.append(S-x[i])

	return out

def Jacobian(a,x,b):
	# First row of our Jacobian
	R = []
	for i in range(5):
		R.append(polyeval(b[i],0))
	out = [R]

	# Creating the other rows
	for i in range(len(x)):

		# Initialising a row
		R = []

		# Loop through the different columns
		for j in range(len(a)):

			# Create the first part of the sum 2*a[j]*p[j](x_i)*p[j]'(x_i)
			S = 2*a[j]*polyeval(b[j],x[i])*polyeval(polydiff(b[j]),x[i])

			for k in range(len(a)):

				if k != j:

					# Add the part a[k]*(p[j](x_i)*p[k]'(x_i)+p[j]'(x_i)*p[k](x_i))
					S += a[k]*(polyeval(polydiff(b[j]),x[i])*polyeval(b[k],x[i])+polyeval(polydiff(b[k]),x[i])*polyeval(b[j],x[i]))

			R.append(S)
		out.append(R)

	return out

# We'll take the different a's to be binary code

def integrate2(f,a,b,solution,N=1000):

	dx = (b-a)/N
	S = f(a,solution)+f(b,solution)

	for i in range(1,N):
		S += 2*f(a+i*dx,solution)


	return 1/2*dx*S

def binary(k):
	out = [0]*5
	for i in range(1,6):

		x = 2**(5-i)
		if k >= x:
			out[-i] = 1
			k += -x

	return out

def e(x,solution):

	return abs(polyeval(solution,x)-math.sqrt(1+x**2))

"""
L = 97
store = 10000
k1 = 0
ap = 0
error = [4,5,6,7,18]
for k in range(L,2**10):
	a = binary(k)
	print("this is k", k)
	if k not in error:
		for i in range(N):

			a = matrix(a)
			J = matrix(Jacobian(a,x,b))
			F = matrix(f(a,x,b))
			
			c = J*a-F

			print(k1)
			a = mpmath.lu_solve(J,c)

		solution = [0]
		for i in range(5):
			solution = polyadd(solution,polysmult(b[i],a[i]))


		s1 = integrate2(e,-1,1,solution)

		if s1 < store:
			store = s1
			ap = a.copy()
			k1 = k
			print(k)
			print(ap)
			print(store)
			
print(store,ap,k1)
"""
a = binary(1)
for i in range(N):

	a = matrix(a)
	J = matrix(Jacobian(a,x,b))
	F = matrix(f(a,x,b))
			
	c = J*a-F

	a = mpmath.lu_solve(J,c)


solution = [0]
for i in range(5):
	solution = polyadd(solution,polysmult(b[i],a[i]))

print(integrate2(e,-1,1,solution))
print(solution)
a = -1
b = 1
points = 100
dx = (b-a)/(points-1)

xlist = [a+i*dx for i in range(points)]
def ef(x):
	return math.sqrt(1+x**2)

eylist = [ef(x) for x in xlist]
plt.plot(xlist,eylist,"-r")

aylist = [polyeval(solution,x) for x in xlist]
plt.plot(xlist,aylist,"-g")
plt.show()

