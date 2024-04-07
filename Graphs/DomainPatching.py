from DESolver5 import *
from datetime import datetime

mpmath.mp.dps = 200
mpmath.mp.prec = 200


N = 25
it = 7
n = 3

D = [0,100]
m = 30
L = (D[1]-D[0])/m

# Initial Conditions
a = 1
b = 0

z1 = []

start = datetime.now()
for k in range(m):

	# Domain of the Patch in the New Space [(k)L,(k+1)L] -> [(k+1)/(k+2),1]
	D1 = [(k+1)/(k+2),1]

	# Initial Guess:
	z = [0,0,3*a-b,b-2*a]

	for i in range(it):
		def g0(x):
			return n*polyeval(z,x)**(n-1)

		def g1(x):
			return 2*x**3/(L**2*(k+1))*(1/(k+1)+1/(x-(k+1)))

		def g2(x):
			return x**4/((k+1)**2*(L)**2)

		def f(x):
			return (n-1)*polyeval(z,x)**n

		g = [g0,g1,g2]

		z = DESolver(g,f,D=D1,XY=[[1],[a]],DXDY=[[1],[b]],basis="c",N=N)

	# Once Complete our List of Solutions gets updated:
	z1.append(z.copy())

	# New Initial Conditions to Ensure Cont. Diff. at the boundary
	a = polyeval(z,(k+1)*L/(L*(k+2)))
	b = (polyeval(polydiff(z),D1[0]))*D1[0]


def h(x):
	k = math.floor(x/L)
	
	return polyeval(z1[k],(k+1)*L/(x+L))

def g(x):
	if x != 0:
		return sin(x/9)/(x/9)

	return 1

def g1(x):
	if x != 0:
		return (sin(x/9)/(x/9))**2

	return 1


points = 300
dx = (D[1]-D[0])/points
xlist = [dx*i+D[0] for i in range(points)]
ylist = [h(x) for x in xlist]
y1list = [g(x) for x in xlist]
y2list = [g1(x) for x in xlist]
plt.plot(xlist,ylist,"-r")
plt.plot(xlist,y1list,"-g")
plt.plot(xlist,y2list,"-b")
plt.show()