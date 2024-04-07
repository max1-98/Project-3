from DESolver5 import *

mpmath.mp.dps = 500
mpmath.mp.prec = 500


N = 25
it = 7
n = 3

D = [0,30]
m = 20
L = (D[1]-D[0])/m

# Initial Conditions
a = 1
b = 0

z1 = []

for k in range(m):

	# Domain of the Patch in the New Space [(k)L,(k+1)L] -> [(k+1)/(k+2),1]
	D1 = [k*L,(k+1)*L]

	# Initial Guess:
	if k != 0:
		z = [(a-b)/(k*L),b]
	else:
		z = [1]

	for i in range(it):
		def g0(x):
			return n*polyeval(z,x)**(n-1)

		def g1(x):
			return 2/x

		def g2(x):
			return 1

		def f(x):
			return (n-1)*polyeval(z,x)**n

		g = [g0,g1,g2]
		z = DESolver(g,f,D=D1,XY=[[k*L],[a]],DXDY=[[k*L],[b]],basis="c",N=N)

	z1.append(z.copy())
	a = polyeval(z,(k+1)*L)
	b = polyeval(polydiff(z),(k+1)*L)

def h(x):
	k = math.floor(x/L)
	
	return polyeval(z1[k],x)

points = 200
dx = (D[1]-D[0])/points
xlist = [dx*i+D[0] for i in range(points)]
ylist = [h(x) for x in xlist]
plt.plot(xlist,ylist,"-r")
plt.show()
