# In this we will solve a simple set of equations, see my Notes for the full equations 


from mpmath import *
from main import * 
mpmath.mp.dps = 100
mpmath.mp.prec = 100

# We'll set this as our number of iterations
N = 50

a = [1,1,0]

x = [mpf('0.92387953251128675612818318939698'), mpf('0.38268343236508977172845998403042')]
b = [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [-1.0, 0.0, 2.0, 0.0, 0.0, 0.0], [1.0, 0.0, -8.0, 0.0, 8.0, 0.0]]

# Creates our f
def f(a,x,b):
	# Each row of f
	out = [a[0]-a[1]+a[2]-1]

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
	out = [[1,-1,1]]

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




for i in range(N):

	a = matrix(a)
	J = matrix(Jacobian(a,x,b))
	F = matrix(f(a,x,b))
	
	c = J*a-F

	a = mpmath.lu_solve(J,c)

print(a)



"""
for i in range(6):

	x = [y[0,0],y[1,0],y[2,0]]
	#J = np.zeros((3,3))

	#JI = np.linalg.inv(J)
	print(y,i)
	y = y - JI.dot(f(x))

"""