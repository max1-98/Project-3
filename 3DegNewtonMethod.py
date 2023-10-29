# In this we will solve a simple set of equations, see my Notes for the full equations 


import numpy as np
# We'll set this as our number of iterations
N = 50

"""
a0 = 1
a1 = 1
a2 = -1
y = np.array([[a0],[a1],[a2]])


def f(x):

	return np.array([[x[0]-1], [-(x[1])**2-2*(x[2])**2-2*x[0]*x[2]+3*x[1]*x[2]+1], [x[1]**2+2*x[2]**2+2*x[0]*x[2]+3*x[1]*x[2]-1]])

for i in range(N):

	x = [y[0,0],y[1,0],y[2,0]]
	J = np.zeros((3,3))

	J[0,0] = 1
	J[1,0] = x[1] - 2*x[2]
	J[1,1] = x[0] - 2*x[1] + 3*x[2]
	J[1,2] = -2*x[0] + 3*x[1] -4*x[2]
	J[2,0] = x[1] + 2*x[2]
	J[2,1] = x[0] + 2*x[1] + 3*x[2]
	J[2,2] = 2*x[0] + 3*x[1] + 4*x[2]
	JI = np.linalg.inv(J)
	print(x,f(x), i)
	print(f(x))
	y = y - JI.dot(f(x))
"""

a0 = 1
a1 = 1
a2 = -1
y = np.array([[a0],[a1],[a2]])


def f(x):

	return np.array([[x[0]-1], [2*(x[1])**2+4*(x[2])**2+2*x[0]*x[1]+4*x[0]*x[2]+6*x[1]*x[2]-1], [1/4*x[1]**2+1/32*x[2]**2+x[0]*x[1]+1/2*x[0]*x[2]+3/16*x[1]*x[2]-1/2]])

for i in range(N):

	x = [y[0,0],y[1,0],y[2,0]]
	J = np.zeros((3,3))

	J[0,0] = 1
	J[1,0] = 2*x[1] - 2*x[2]
	J[1,1] = 2*x[0] + 4*x[1] + 6*x[2]
	J[1,2] = 4*x[0] + 6*x[1] + 8*x[2]
	J[2,0] = x[1] + 1/2*x[2]
	J[2,1] = x[0] + 1/2*x[1] + 3/16*x[2]
	J[2,2] = 1/2*x[0] + 3/16*x[1] + 1/16*x[2]
	JI = np.linalg.inv(J)
	print(x,f(x), i)
	print(f(x))
	y = y - JI.dot(f(x))

