from DESolver5 import *

mpmath.mp.dps = 300
mpmath.mp.prec = 300


plt.clf()

N = 25
it = 35
D = [0,25]


n = 5
points = 200
dx = (D[1]-D[0])/(points-1)
xlist = [D[0]+i*dx for i in range(points)]
y1list = [1/sqrt(1+x**2/3) for x in xlist]


z1 = [1]
for i in range(it):

	# OUr ODE's Coefficients
	def g0(x):
		return n*polyeval(z1,x)**(n-1)

	def g1(x):
		return 2/x

	def g2(x):
		return 1

	# What our ODE equals
	def f(x):
		return (n-1)*polyeval(z1,x)**n

	g = [g0,g1,g2]
	z2 = DESolver(g,f,D=D,XY=[[0],[1]],DXDY=[[0],[0]],basis="l",N=N)

	def ser(x):
		return abs(x*n*polyeval(z1,x)**(n-1)*polyeval(z2,x)+2*polyeval(polydiff(z2),x)+x*polyeval(polydiffn(z2,2),x)-(n-1)*polyeval(z1,x)**n)


	z1 = z2

	ay1list = [polyeval(z1,x) for x in xlist]
	lab = "Iteration: " + str(i+1)

	plt.plot(xlist,ay1list,":g",label="Iterate")
	plt.plot(xlist,y1list,"-r",label="Exact Solution")
	plt.legend(loc="upper left")
	ax = plt.gca()
	plt.xlabel("x")
	plt.ylabel("y")
	plt.title("Iteration "+str(i+1))

	ax.set_ylim([-0.2, 1.2])
	plt.savefig(" IT "+ str(i+1))
	plt.clf()