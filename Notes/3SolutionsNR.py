import math
import matplotlib.pyplot as plt
from main import *
import mpmath
import latex

mpmath.mp.dps = 100
mpmath.mp.prec = 100


def g(x):
	return math.sqrt(1+x**2)


plt.grid()


p = [[ 0.50345836984584010841367435832,-0.88889486900965212207079693309,-0.39235323885549223048447129142],[ -1.289403497522859017886750519,-1.7174826862356650298244754283,0.57192081128719398806227509065],[    1.2173411583055409963875616632, 0.20842027928171569810339479336, -0.0089208790238252982841668698426]]

plt.xlabel(r"$x$")
plt.ylabel(r'$y$')

B = [[1],[-1,0,2],[1,0,-8,0,8]]

def f(x,p,B,i):

	S = 0
	N = 0
	for j in range(len(B)):
		S += polyeval(B[j],x)*p[i][j]
		#N += polyeval(polydiff(B[j]),x)*p[i][j]

	return S
	return S*N

a = -1
b = 1
points = 100
dx = (b-a)/points
xlist = [-1+dx*i for i in range(points)]
ylist = []
for i in range(3):
	ylist1 = []

	for x in xlist:
		ylist1.append(f(x,p,B,i))

	ylist.append(ylist1)

c = ["b","y","#4CAF50"]

aylist = []
for x in xlist:
	aylist.append(g(x))


plt.plot(xlist,aylist, "-r",label=r"$y=\sqrt{1+x^2}$")
#plt.plot(xlist,xlist,"-", c="r",label=r"$y=x$")
plt.plot(xlist,ylist[0],"-"+c[0],label=r"$y=y_1$")
plt.plot(xlist,ylist[1],"-"+c[1],label=r"$y=y_2$")
plt.plot(xlist,ylist[2],"--", c=c[2],label=r"$y=y_3$")
plt.legend(loc="upper left")


# Two Collocation points
#plt.plot(0.923879,0.923879,"o")
#plt.plot(0.382683,0.382683,"o")

# Initial Point
plt.plot(0,1,"o")

def g(x,p,B,i):

	return abs(f(x,p,B,i)-x)


for i in range(3):
	

	N = 1000
	dx = (b-a)/N
	S1 = g(a,p,B,i)+g(b,p,B,i)

	for j in range(1,N):
		S1 += 2*g(a+j*dx,p,B,i)


	print(1/2*dx*S1)
	# Integrate abs(f(x,p,B,i)-x)


plt.show()