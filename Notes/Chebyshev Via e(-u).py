from DESolver5 import *

mpmath.mp.dps = 500
mpmath.mp.prec = 500


plt.clf()
plt.grid()
plt.title("Laguerre Iterates vs Exact Solution on [0,30]")
N = 50
it = 100



D = [exp(-30),1]

n = 3
points = 200
dx = (D[1]-D[0])/(points-1)
xlist = [D[0]+i*dx for i in range(points)]

plt.clf()
plt.grid()

z1 = [3/2,-1,1/2]
for i in range(it):

	# OUr ODE's Coefficients
	def g0(x):
		return n*polyeval(z1,x)**(n-1)

	def g1(x):
		return (x*log(x)+2*x)/log(x)

	def g2(x):
		return x**2

	# What our ODE equals
	def f(x):
		return (n-1)*polyeval(z1,x)**n

	g = [g0,g1,g2]
	z2 = DESolver(g,f,D=D,XY=[[1],[1]],DXDY=[[1],[0]],basis="c",N=N)
	z1 = z2

	ay1list = [polyeval(z1,x) for x in xlist]
	lab = "Iteration: " + str(i+1)
	if i < it-1:
		if i == 0:
			plt.plot(xlist,ay1list,":g",label="Iterates 1 to "+str(it-1))
		else:
			plt.plot(xlist,ay1list,":g")
	else:
		plt.plot(xlist,ay1list,"-r",label="Final Solution")

print(z2)
plt.legend()
ax = plt.gca()

ax.set_ylim([-10, 10])
plt.show()


