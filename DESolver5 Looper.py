from DESolver5 import *

mpmath.mp.dps = 1000
mpmath.mp.prec = 1000

#DESolver(g,f,D,XY=[],DXDY=[], basis="c", N=20)




z = [0,1]
D = [-1,1]


points = 200
dx = (D[1]-D[0])/(points-1)
xlist = [D[0]+i*dx for i in range(points)]

def h(x):
	return math.asin(x)

eylist = [h(x) for x in xlist]

plt.plot(xlist,eylist,"-k",label="Exact Solution", alpha=0.8)

for i in range(7):

	def g0(x):
		return -math.sin(polyeval(z,x))*polyeval(polydiff(z.copy()),x)

	def g1(x):
		return math.cos(polyeval(z,x))

	

	# What our ODE equals
	def f(x):

		return 1-math.sin(polyeval(z,x))*polyeval(polydiff(z.copy()),x)*polyeval(z,x)

	g = [g0,g1]
	yi1 = z.copy()
	z = DESolver(g,f,D=D,XY=[[0],[0]],DXDY=[],basis="c",N=50)
	yi = z.copy()



	#print("This is solution #", i, ": ",z)
	aylist = [polyeval(z,x) for x in xlist]
	lab = "Iteration: "+ str(i+1)
	if i == 3:
		p1 = z.copy()

	if i == 5:
		p2 = z.copy()

	if i%4==1:
		plt.plot(xlist,aylist,":",label=lab)

plt.title("Domain $[-1,1]$")
plt.grid()
plt.legend()
plt.show()
plt.clf()

a = 0.9
b = 1
dx = (1-0.9)/(points-1)
xlist = [0.9+i*dx for i in range(points)]
y1list = [polyeval(p1,x) for x in xlist]
y2list = [polyeval(p2,x) for x in xlist]
eylist = [h(x) for x in xlist]

plt.plot(xlist,eylist,"-k",label="Exact Solution")
plt.plot(xlist, y1list, "-r", label="Iteration: 2")
plt.plot(xlist,y2list,"-g", label="Iteration: 6")


plt.title("Domain $[0.9,1]$")
plt.grid()
plt.legend()

plt.show()





