


# In this code we will use extremely basic Spectral Methods to approximate the solution
# to the Simple Harmonic Oscillator (d^2/dt^2)x+x=0 x(0)=1, x(pi/2)=0 this has exact solution x(t)=cost)


from math import pi
#import matplotlib.pyplot
#import numpy


def f(t,a):

	s = 0
	for i in range(len(a)):
		s += a[i]*t**i

	return s

# N is the number of coefficients

# Here we are forming the points which we want to be exact
N=2
dx = (pi/2-1)/N
xl = [dx*i for i in range(N)]

# Next we create our matrix that we need to solve

M = []

for i in range(N):
	R = [1, xl[i]]
	for j in range(2,N):
		R.append(j*(j-1)*(xl[i])**(j-2)+xl[i]**j)
	M.append(R)


print(M)




