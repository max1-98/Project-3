

# in this we are going to solve y''+y=x^2 with y(0)=-1 and y(pi/2)=(pi/2)^2-2


import numpy as np
import math
import matplotlib.pyplot as plt



# This is our exact solution 
def ef(x):
	return math.cos(x)+x**2-2


points = 100
dx = (math.pi/2)/points
xlist = [dx*i for i in range(points+1)]
ylist = [ef(x) for x in xlist]
plt.plot(xlist,ylist, "-r")


# This will be the number of orthogonal polynomials we'll use

N = 10

M = []
# X values of the initial Conditions
X = [0, math.pi/2]
Y = [-1, (math.pi/2)**2-2]
for i in range(2):
	R=[]
	for j in range(N):
		R.append(X[i]**j)
	M.append(R)


DX = (X[1]-X[0])/(N-1)
x = [DX*(i+1) for i in range(N-2)]
print(x)

for i in range(N-2):
	R = [1,x[i]]
	for j in range(2,N):
		R.append(j*(j-1)*x[i-2]**(j-2)+x[i-2]**j)

	M.append(R)

B = Y.copy()

for i in range(N-2):
	B.append(x[i]**2)

M = np.array(M)
B = np.array(B)
print(M)

a = np.linalg.solve(M,B)
print(a)

def af(x,a):
	S = 0
	for i in range(len(a)):
		S += a[i]*x**i

	return S
print(B)

yalist = [ af(x,a) for x in xlist]
plt.plot(xlist,yalist, "-g")
plt.show()