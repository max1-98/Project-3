import math
import matplotlib.pyplot as plt
from main import *
import mpmath
mpmath.mp.dps = 100
mpmath.mp.prec = 100


"""
errors = [0.0048917617858682155724726426324210895164765038248357351712768154136054779628087618673951221697629128176394142407013930082147581611187713536763031233168796969624906272258112149192473260348091253241207299652499150226011968687031056523182762726928450136894103570595317012154597129778525657781628691728182,0.000064299471263961597456263153978,0.00000096762504846183790060846744617,0.00000096603337382751210508769867476,0.00000096603342022216059782360628558,0.00000096603342022216059790936142282]


e = []

for error in errors:
	e.append(-math.log(error,10))

plt.plot(e,"-or")

plt.show()
"""

p = [[ 0.50345836984584010841367435832,-0.88889486900965212207079693309,-0.39235323885549223048447129142],[ -1.289403497522859017886750519,-1.7174826862356650298244754283,0.57192081128719398806227509065],[    1.2173411583055409963875616632, 0.20842027928171569810339479336, -0.0089208790238252982841668698426]]
B = [[1],[-1,0,2],[1,0,-8,0,8]]

def f(x,p,B,i):

	S = 0
	N = 0
	for j in range(len(B)):
		S += polyeval(B[j],x)*p[i][j]
		N += polyeval(polydiff(B[j]),x)*p[i][j]

	#return S
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

c = ["b","y","g"]

for i in range(3):
	plt.plot(xlist,ylist[i],"-"+c[i])

"""
def g(x):
	return math.sqrt(1+x**2)

aylist = []
for x in xlist:
	aylist.append(g(x))

plt.plot(xlist,aylist, "-r")
"""
plt.plot(xlist,xlist,"-m")

plt.plot(0.923879,0.923879,"o")
plt.plot(0.382683,0.382683,"o")
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

