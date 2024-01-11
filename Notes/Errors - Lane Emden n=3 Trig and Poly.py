from DESolver5 import *

mpmath.mp.dps = 1000
mpmath.mp.prec = 1000


plt.clf()
plt.grid()

N = 63
c = [1]
s = []


z = [c,s]
D = [0,3]

n = 3
points = 200
dx = (D[1]-D[0])/(points-1)
xlist = [D[0]+i*dx for i in range(points)]


for i in range(6):

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

# Best Solution Found
p2 = [mpf('1.0000000000000000000000000000032'), mpf('1.6543738222286355977750850526805e-31'), mpf('-0.16666417267038417201429394558274'), mpf('-0.000082957743306267661150973724981338'), mpf('0.025605131576298947177404546425636'), mpf('-0.0020693138682529209065091112245346'), mpf('0.00033744217231567327571755035623571'), mpf('-0.0051820675069185669710686712039039'), mpf('0.0049002297718746633097810693107764'), mpf('-0.0024172645012616679802277913357769'), mpf('0.00078588078352131481990817087547127'), mpf('-0.0001829408346829232838159991255107'), mpf('0.000031604756017320192649462133458241'), mpf('-0.0000041081849558739600674186262096809'), mpf('0.00000040154867323610561326153897040344'), mpf('-0.000000029115664645186339119075126205264'), mpf('0.0000000015204497841228067678574927410511'), mpf('-5.4095346916174923801327285037428e-11'), mpf('1.1742743039932198914195953196549e-12'), mpf('-1.1737966146828494002811853010005e-14')]

z1 = [1]
for i in range(6):

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

	print(-math.log(integrate(ser,0,D[1]),10)/D[1], "Polynomial Error, Iteration: ", i)
	z1 = z2

	ay1list = [polyeval(z1,x) for x in xlist]
	lab = "Iteration: " + str(i+1)


# Exact Error 
def eer(x):
	return abs(polyeval(z1,x)-trigeval(z,x))
# Residual Error
def rer(x):
	return abs(x*trigeval(z,x)**(n)+2*trigeval(trigdiff(z),x)+x*trigeval(trigdiffn(z,2),x))

def rerp(x):
	return abs(x*polyeval(z1,x)**(n)+2*polyeval(polydiff(z1),x)+x*polyeval(polydiffn(z1,2),x))

print(-math.log(integrate(eer,0,2),10)/D[1], " Exact Error")
print(-math.log(integrate(rer,0,2),10)/D[1], " Residual Error Trig")
print(-math.log(integrate(rerp,0,2),10)/D[1], " Residual Error Polynomial")
plt.plot(xlist,aylist,":",label=lab)
eylist = [polyeval(p2,x) for x in xlist]
plt.plot(xlist,eylist,"-g")
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