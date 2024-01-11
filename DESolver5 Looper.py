from DESolver5 import *

mpmath.mp.dps = 500
mpmath.mp.prec = 500


plt.clf()
plt.grid()

N = 25
it = 7
c = [1]
s = []


z = [c,s]
D = [0,10]

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


	print(-math.log(integrate(ser,0,D[1]),10)/D[1], "Laguerre Error, Iteration: ", i)
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

plt.legend()
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
plt.legend()
plt.show()


###
"""

N = 13
-0.23298537719997933 Trigonometric Error, Iteration:  0
-0.09550281372765011 Trigonometric Error, Iteration:  1
-0.062414257897812554 Trigonometric Error, Iteration:  2
-0.061971073732746455 Trigonometric Error, Iteration:  3
-0.06197102790718715 Trigonometric Error, Iteration:  4
-0.061971027907187 Trigonometric Error, Iteration:  5
-0.23299011753313312 Polynomial Error, Iteration:  0
-0.09547042148808955 Polynomial Error, Iteration:  1
-0.06235174549876707 Polynomial Error, Iteration:  2
-0.061907670066403186 Polynomial Error, Iteration:  3
-0.0619076225273964 Polynomial Error, Iteration:  4
-0.06190762252739613 Polynomial Error, Iteration:  5
1.428290865700191  Exact Error
0.8324844158452686  Residual Error Trig
1.9206957849312234  Residual Error Polynomial

N = 23
-0.23299011625728783 Trigonometric Error, Iteration:  0
-0.09547040922486921 Trigonometric Error, Iteration:  1
-0.06235175150349076 Trigonometric Error, Iteration:  2
-0.06190767603328138 Trigonometric Error, Iteration:  3
-0.06190762849492646 Trigonometric Error, Iteration:  4
-0.061907628494926197 Trigonometric Error, Iteration:  5
-0.23299011725715502 Polynomial Error, Iteration:  0
-0.0954704012312339 Polynomial Error, Iteration:  1
-0.062351736129308506 Polynomial Error, Iteration:  2
-0.061907660549862825 Polynomial Error, Iteration:  3
2.7126256026953555  Exact Error
1.9128199388065437  Residual Error Trig
3.5971059861654475  Residual Error Polynomial

-0.19109780022903847 Polynomial Error, Iteration:  0
-0.14024298292542375 Polynomial Error, Iteration:  1
-0.09180318961133646 Polynomial Error, Iteration:  2
-0.09421768395046881 Polynomial Error, Iteration:  3
-0.06156375604748012 Polynomial Error, Iteration:  4
-0.04327087545189924 Polynomial Error, Iteration:  5
-0.03931008259315026 Polynomial Error, Iteration:  6

-0.19138138525594317 Polynomial Error, Iteration:  0
-0.13964281919652366 Polynomial Error, Iteration:  1
-0.09179464094830983 Polynomial Error, Iteration:  2
-0.07874667929518216 Polynomial Error, Iteration:  3
-0.03379899920565697 Polynomial Error, Iteration:  4
-0.02268107584645695 Polynomial Error, Iteration:  5
-0.0227573849768559 Polynomial Error, Iteration:  6

"""

"""
N = 33
-0.23299011725699745 Trigonometric Error, Iteration:  0
-0.09547040123254436 Trigonometric Error, Iteration:  1
-0.06235173613180109 Trigonometric Error, Iteration:  2
-0.061907660552370535 Trigonometric Error, Iteration:  3
-0.06190761301384781 Trigonometric Error, Iteration:  4
-0.061907613013847514 Trigonometric Error, Iteration:  5
-0.23299011725715502 Polynomial Error, Iteration:  0
-0.09547040123093847 Polynomial Error, Iteration:  1
-0.06235173613249722 Polynomial Error, Iteration:  2
-0.06190766055304588 Polynomial Error, Iteration:  3
-0.06190761301452252 Polynomial Error, Iteration:  4
-0.06190761301452225 Polynomial Error, Iteration:  5
3.984074126903842  Exact Error
3.064936557440154  Residual Error Trig
4.698440412769904  Residual Error Polynomial


N = 43 
-0.2329901172571551 Trigonometric Error, Iteration:  0
-0.09547040123120425 Trigonometric Error, Iteration:  1
-0.06235173612921591 Trigonometric Error, Iteration:  2
-0.0619076605497691 Trigonometric Error, Iteration:  3
-0.06190761301124639 Trigonometric Error, Iteration:  4
-0.061907613011246095 Trigonometric Error, Iteration:  5
-0.23299011725715502 Polynomial Error, Iteration:  0
-0.09547040117958833 Polynomial Error, Iteration:  1
-0.06235173615308683 Polynomial Error, Iteration:  2
-0.06190766057677016 Polynomial Error, Iteration:  3
-0.06190761303826454 Polynomial Error, Iteration:  4
-0.06190761303826426 Polynomial Error, Iteration:  5
5.232294516113281  Exact Error
4.249255467960627  Residual Error Trig
4.763390804631022  Residual Error Polynomial

N = 53
-0.23299011725715513 Trigonometric Error, Iteration:  0
-0.09547040123120408 Trigonometric Error, Iteration:  1
-0.06235173612921561 Trigonometric Error, Iteration:  2
-0.06190766054976881 Trigonometric Error, Iteration:  3
-0.061907613011246095 Trigonometric Error, Iteration:  4
-0.0619076130112458 Trigonometric Error, Iteration:  5
-0.23299011725715502 Polynomial Error, Iteration:  0
-0.09547040106707193 Polynomial Error, Iteration:  1
-0.062351740706108306 Polynomial Error, Iteration:  2
-0.061907664892301824 Polynomial Error, Iteration:  3
-0.06190761735277303 Polynomial Error, Iteration:  4
-0.06190761735277275 Polynomial Error, Iteration:  5
5.397767873960191  Exact Error
5.005569156510812  Residual Error Trig
4.551781033056177  Residual Error Polynomial

N = 63
-0.23299011725715513 Trigonometric Error, Iteration:  0
-0.09547040123120408 Trigonometric Error, Iteration:  1
-0.06235173612921561 Trigonometric Error, Iteration:  2
-0.06190766054976881 Trigonometric Error, Iteration:  3
-0.061907613011246095 Trigonometric Error, Iteration:  4
-0.0619076130112458 Trigonometric Error, Iteration:  5
-0.23299011725715502 Polynomial Error, Iteration:  0
-0.09547040065414131 Polynomial Error, Iteration:  1
-0.06235141678596786 Polynomial Error, Iteration:  2
-0.061907386808409036 Polynomial Error, Iteration:  3
-0.06190733864709338 Polynomial Error, Iteration:  4
-0.06190733864580908 Polynomial Error, Iteration:  5
5.27490328914896  Exact Error
4.964472935654925  Residual Error Trig
4.29703468584634  Residual Error Polynomial


"""