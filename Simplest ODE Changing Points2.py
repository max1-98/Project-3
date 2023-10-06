# Solving dy/dx = -y with y(0) = 1, but in this we set N=3 and change the initial points x_0 and x_1 to try optimise the rate of convergence.




import math
import numpy as np
import matplotlib.pyplot as plt

# Our exact function and our approximate below for coefficients a
def ef(x):

	return math.e**(-x)

def af(x,a):

	S = 0

	for i in range(len(a)):
		S += a[i]*x**i

	return S


# This function is the solution to our GJ matrix
def coefficients(x):
	t = x[0]
	u = x[1]

	return [1, ((u)*(u-2)+t*(2-t))/(u*(u-2)*(t-1)-t*(t-2)*(u-1)),(u-t)/(u*(u-2)*(t-1)-t*(t-2)*(u-1))]

# We will start this off like a plot but really we will be using these points to form an approximation for the integral of y_N-e^x

points = 100
dx = (2/(points-1))
xlist = [-1+i*dx for i in range(points)]
eylist = [ef(x) for x in xlist]

# Initial the error to be massive
store = 9999

# We'll let the difference begin as the same as dx
DX = dx

# Likewise original a and b will be -1 and -1
b = -1
c = -1


# The i and j are so that we loop through all the different 
for z in range(3):
	for i in range(points):
		for j in range(points):

			# Makes sure two of the initial points aren't equal otherwise we will get no sol in our GJ matrix
			if i != j:

				# Creating our two x values

				x = [b+DX*i, c+DX*j]



				a = coefficients(x)

				# Trapezium Rule to find the error
				S = abs(af(-1,a)-eylist[0])+abs(af(1,a)-eylist[-1])

				for k in range(1,len(xlist)-1):

					S += abs(2*(af(xlist[k],a)-eylist[k]))


				if S < store:
					store = S
					x1 = x[0]
					x2 = x[1]
	print("This is cycle number: ", z)
	print("Lowest error:", 1/2*dx*store)
	print("The x values producing this: ", x1, " and ", x2)

	DX = dx/50
	b = x1
	c = x2



print("Lowest error:", 1/2*dx*store)
print("The x values producing this: ", x1, " and ", x2)
"""



def coefficients1(x):
	N=3
	# Initialising our matrices in Ma=B, R is going to be our row of M
	M = []
	B = [1]
	R = [1]

	# Creating B and the first row of M
	for i in range(N-1):
		R.append(0)
		B.append(0)

	M.append(R)



	for n in range(2):

		R = []
		for k in range(3):
			R.append(k*x[n]**(k-1)+x[n]**k)

		
		M.append(R)


	M = np.array(M)
	B = np.array(B)

	a = np.linalg.solve(M,B)

	return a


print(coefficients([-1,1]))
print(coefficients1([-1,1]))
"""