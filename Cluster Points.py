
# Solving dy/dx = -y with y(0) = 1

import math
import numpy as np
import matplotlib.pyplot as plt

def ef(x):

	return math.e**(-x)
points = 100
dx = 2/(points-1)
xlist = [-1+i*dx for i in range(points)]
eylist = [ef(x) for x in xlist]

plt.plot(xlist,eylist,"-r")


N = 3

# Our Matrix in Mx=B where
M = []
B = [1]
R = [1]

for i in range(N-1):
	R.append(0)
	B.append(0)

M.append(R)

x = [-0.3668341708542714,0.44723618090452266]


for i in range(N-1):

	R = []
	for j in range(N):
		R.append(j*x[i]**(j-1)+x[i]**j)

	M.append(R)

M = np.array(M)
B = np.array(B)

a = np.linalg.solve(M,B)

print(a)
def af(x,a):

	S = 0

	for i in range(len(a)):
		S += a[i]*x**i

	return S


a1ylist = [af(x,a) for x in xlist]


#plt.plot(xlist,a1ylist,"-g")
plt.plot(xlist,a1ylist,"-b")


plt.show()