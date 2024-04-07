from DESolver5 import *

mpmath.mp.dps = 300
mpmath.mp.prec = 300


plt.clf()

N = 25
it = 5
D = [0,4]
D1 = [1/(1+D[1]),1]


n = 5
points = 200
dx = (D[1]-D[0])/(points-1)
xlist = [D[0]+i*dx for i in range(points)]
y1list = [1/sqrt(1+x**2/3) for x in xlist]

plt.plot(xlist,y1list,"-r",label=r"$y=\frac{1}{\sqrt{1+x^2/3}}$")
z1 = [0,2,-1]
for i in range(it):

	# Our ODE's Coefficients
	def g0(x):
		return n*polyeval(z1,x)**(n-1)

	def g1(x):
		return 2*x**4/(x-1)

	def g2(x):
		return x**4

	# What our ODE equals
	def f(x):
		return (n-1)*polyeval(z1,x)**n


	g = [g0,g1,g2]
	z2 = DESolver(g,f,D=D,XY=[[1],[1]],DXDY=[[1],[0]],basis="c",N=N)

	z1 = z2.copy()

	def h(x):
		return polyeval(z1,1/(1+x))

	# Residual Error
	def err(x):
		return abs((x-1)*x**4*polyeval(polydiffn(z1,2),x)+2*x**4*polyeval(polydiff(z1),x)+(x-1)*polyeval(z1,x)**n)*1/(x)**2

	print(-math.log(1/(D[1]-D[0])*integrate(err,D1[0],D1[1]),10), " IT" + str(i+1))

	ay1list = [h(x) for x in xlist]

	
	
	
	plt.plot(xlist,ay1list,"--",label=r"$y_{"+str(i+1)+"}$")
	plt.legend(loc="upper right")

	# plt.plot(xlist,y1list,"-r",label="Exact Solution")
	#plt.plot(xlist,ay1list,":g",label="Iterate")
	#plt.title("Iteration "+str(i+1))
	#plt.savefig(" IT "+ str(i+1))
	#plt.clf()

ax = plt.gca()
plt.xlabel("x")
plt.ylabel("y")
plt.title(r"Rational Chebyshev Approx. Sol. vs Exact Sol. (Lane Emden $n=5$ on $[0,4]$)")
plt.plot(xlist,y1list,"-r",label=r"$y=\frac{1}{\sqrt{1+x^2/3}}")
ax.set_ylim([0.3, 1.1])
plt.show()