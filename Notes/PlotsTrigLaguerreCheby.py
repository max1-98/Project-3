from DESolver5 import *

mpmath.mp.dps = 500
mpmath.mp.prec = 500


plt.clf()
plt.grid()

N = 30
it = 7
c = [1]
s = []


z = [c,s]
D = [0,30]

n = 3
points = 200
dx = (D[1]-D[0])/(points-1)
xlist = [D[0]+i*dx for i in range(points)]

""" Trig Functions
for i in range(it):

	# OUr ODE's Coefficients
	def g0(x):
		return n*trigeval(z,x)**(n-1)

	def g1(x):
		return 2/x

	def g2(x):
		return 1

	# What our ODE equals
	def f(x):
		return (n-1)*trigeval(z,x)**n

	g = [g0,g1,g2] 
	yi1 = z.copy()
	z2 = DESolver(g,f,D=D,XY=[[0],[1]],DXDY=[[0],[0]],basis="t",N=N)
	yi = z.copy()

	def ser(x):
		return abs(x*n*trigeval(z,x)**(n-1)*trigeval(z2,x)+2*trigeval(trigdiff(z2),x)+x*trigeval(trigdiffn(z2,2),x)-(n-1)*trigeval(z,x)**n)
	print(-math.log(integrate(ser,0,D[1]),10)/D[1], "Trigonometric Error, Iteration: ", i)

	z = z2

	aylist = [trigeval(z,x) for x in xlist]
	lab = "Iteration: "+ str(i+1)
	plt.plot(xlist,aylist,":r",label=lab)
"""
plt.clf()
plt.grid()
plt.title("Laguerre Basis")

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
	yi1 = z.copy()
	z2 = DESolver(g,f,D=D,XY=[[0],[1]],DXDY=[[0],[0]],basis="l",N=N)
	yi = z.copy()

	def ser(x):
		return abs(x*n*polyeval(z1,x)**(n-1)*polyeval(z2,x)+2*polyeval(polydiff(z2),x)+x*polyeval(polydiffn(z2,2),x)-(n-1)*polyeval(z1,x)**n)


	#print(-math.log(integrate(ser,0,D[1]),10)/D[1], "Laguerre Error, Iteration: ", i)
	z1 = z2

	ay1list = [polyeval(z1,x) for x in xlist]
	lab = "Iteration: " + str(i+1)
	if i < it:
		if i == 0:

			plt.plot(xlist,ay1list,":g",label="Iterates")
		else:
			plt.plot(xlist,ay1list,":g")

	if i == it-1:
		plt.plot(xlist,ay1list,"-r",label="Final Approximate Solution")
print(z1)
plt.legend()
ax = plt.gca()

ax.set_ylim([-10, 10])
plt.show()


plt.clf()

plt.title("Chebyshev Basis")
plt.grid()

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
	yi1 = z.copy()
	z2 = DESolver(g,f,D=D,XY=[[0],[1]],DXDY=[[0],[0]],basis="c",N=N)
	yi = z.copy()

	def ser(x):
		return abs(x*n*polyeval(z1,x)**(n-1)*polyeval(z2,x)+2*polyeval(polydiff(z2),x)+x*polyeval(polydiffn(z2,2),x)-(n-1)*polyeval(z1,x)**n)


	print(-math.log(integrate(ser,0,D[1]),10)/D[1], "Chebyshev Error, Iteration: ", i)
	z1 = z2

	ay1list = [polyeval(z1,x) for x in xlist]
	lab = "Iteration: " + str(i+1)
	if i < it:
		if i ==0:
			plt.plot(xlist,ay1list,":y",label="Iterates")
		else:
			plt.plot(xlist,ay1list,":y")

	if i == it-1:
		plt.plot(xlist,ay1list,"-r",label="Final Approximate Solution")
print(z1)
plt.legend()
plt.show()
