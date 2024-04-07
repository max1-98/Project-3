from DESolver5 import *

mpmath.mp.dps = 100
mpmath.mp.prec = 100

#DESolver(g,f,D,XY=[],DXDY=[], basis="c", N=20)
M=3


plt.title(r"Chebyshev approx. sol. $y=y_i$ (Lane-Emden $n=3$ on $[0,10]$)")
z = [1]
D = [0,10]


points = 100
dx = (D[1]-D[0])/(points-1)
xlist = [D[0]+i*dx for i in range(points)]

plt.plot(xlist,[1 for x in xlist],"--",label=r"$y=y_0$")

for i in range(20):

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
	z = DESolver(g,f,D=D,XY=[[0],[1]],DXDY=[[0],[0]],basis="c",N=25)
	yi = z.copy()

	def h(x):

		return abs(x*polyeval(polydiff(polydiff(yi.copy())),x)+2*polyeval(polydiff(yi.copy()),x)+x*polyeval(yi.copy(),x)**M)

	a = integrate(h,0,D[1],N=1000)

	print(-math.log(1/(D[1]-D[0])*a,10),"Residual Error - ODE")


	aylist = [polyeval(z,x) for x in xlist]
	plt.plot(xlist,aylist,"--",label=r"$y=y_{"+str(i+1)+"}$")

plt.legend()
plt.show()
plt.clf()