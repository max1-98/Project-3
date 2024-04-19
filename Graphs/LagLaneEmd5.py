from DESolver5 import *

mpmath.mp.dps = 300
mpmath.mp.prec = 300

#DESolver(g,f,D,XY=[],DXDY=[], basis="c", N=20)
M=5


plt.title(r"Laguerre approx. sol. $y=y_{50}$, vs. Exact sol. $y=\frac{1}{\sqrt{1+x^2/3}}$. (Lane-Emden $n=5$ on $[0,30]$)")
z = [1]
D = [0,30]


points = 100
dx = (D[1]-D[0])/(points-1)
xlist = [D[0]+i*dx for i in range(points)]

plt.plot(xlist,[1/sqrt(1+x**2/3) for x in xlist],"-r",label=r"$y=\frac{1}{\sqrt{1+x^2/3}}$")
#plt.plot(xlist,[1 for x in xlist],"--",label=r"$y=y_0$")

for i in range(50):

	def g0(x):
		return M*polyeval(z,x)**(M-1)

	def g1(x):
		return 2/x

	def g2(x):
		return 1

	# What our ODE equals
	def f(x):

		return (M-1)*polyeval(z,x)**(M)

	g = [g0,g1,g2]
	yi1 = z.copy()
	z = DESolver(g,f,D=D,XY=[[0],[1]],DXDY=[[0],[0]],basis="l",N=25)
	yi = z.copy()

	aylist = [polyeval(z,x) for x in xlist]
	#plt.plot(xlist,aylist,"--",label=r"$y=y_{"+str(i+1)+"}$")

ax = plt.gca()
plt.xlabel("x")
plt.ylabel("y")
ax.set_ylim([-0.2, 1.2])
plt.plot(xlist,aylist,"--",label=r"$y=y_{"+str(i+1)+"}$")
plt.legend()
plt.show()
plt.clf()