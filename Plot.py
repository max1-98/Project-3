import math
import matplotlib.pyplot as plt
from main import *
import mpmath
import latex
from mpmath import *

mpmath.mp.dps = 100
mpmath.mp.prec = 100


def g(x):
	return math.sqrt(1+x**2)

def h(x,s):

	return polyeval(s,x)

plt.title("On the intveral " + r"$[0.999,1]$")
plt.grid()
plt.xlabel(r"$x$")
plt.ylabel(r'$y$')

p = [[ 0.50345836984584010841367435832,-0.88889486900965212207079693309,-0.39235323885549223048447129142],[ -1.289403497522859017886750519,-1.7174826862356650298244754283,0.57192081128719398806227509065],[    1.2173411583055409963875616632, 0.20842027928171569810339479336, -0.0089208790238252982841668698426]]



B = [[1],[-1,0,2],[1,0,-8,0,8]]

def f(x,p,B,i):

	S = 0

	for j in range(len(B)):
		S += polyeval(B[j],x)*p[i][j]


	return S


a = 0.999
b = 1
points = 200
dx = (b-a)/points
xlist = [a+dx*i for i in range(points)]
ylist = []

for x in xlist:
	ylist.append(f(x,p,B,2))

s = [mpf('0.99999999999999999999999999999921'), mpf('0.0'), mpf('0.49973666857582575077137573566055'), mpf('0.0'), mpf('-0.12062258563421877376961495106674'), mpf('0.0'), mpf('0.045535878902564385505941386762077'), mpf('0.0'), mpf('-0.010425003185446731647617695197307'), mpf('0.0')]
ylist1 = [h(x,s) for x in xlist]

c = ["b","#4CAF50","y"]

aylist = []
for x in xlist:
	aylist.append(g(x))


plt.plot(xlist,aylist, "-",c="#000",label=r"$y=\sqrt{1+x^2}$")
plt.plot(xlist,ylist,":", c=c[0],label=r"$y=y_3$")
plt.plot(xlist,ylist1,":",c="r",label=r"$y=y_{10}$")
plt.legend(loc="upper left")

plt.show()