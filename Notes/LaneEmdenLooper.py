from DESolver5 import *

mpmath.mp.dps = 100
mpmath.mp.prec = 100

#DESolver(g,f,D,XY=[],DXDY=[], basis="c", N=20)
M=3



z = [1]
D = [0,1/3]


points = 100
dx = (D[1]-D[0])/(points-1)
xlist = [D[0]+i*dx for i in range(points)]



for i in range(5):

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
	z = DESolver(g,f,D=D,XY=[[0],[1]],DXDY=[[0],[0]],basis="c",N=40)
	yi = z.copy()

	def h(x):

		return abs(x*polyeval(polydiff(polydiff(yi.copy())),x)+2*polyeval(polydiff(yi.copy()),x)+x*polyeval(yi.copy(),x)**3)

	a = integrate(h,0,1/3,N=1000)

	print(-math.log(6*a,10),"Differential Error - ODE")

	def k(x):
		return abs(x*polyeval(polydiff(polydiff(yi.copy())),x)+2*polyeval(polydiff(yi.copy()),x)+3*x*(polyeval(yi1.copy(),x)**2)*polyeval(yi.copy(),x)-2*polyeval(yi1,x)**3)

	b = integrate(k,0,1/3,N=1000) 
	print(-math.log(6*b,10),"Differential Error - NKODE")

	#print("This is solution #", i, ": ",z)
	aylist = [polyeval(z,x) for x in xlist]
	lab = "Iteration: "+ str(i)
	plt.plot(xlist,aylist,"--",label=lab)
	
plt.show()
plt.clf()