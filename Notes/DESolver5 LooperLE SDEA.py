from DESolver5 import *

mpmath.mp.dps = 100
mpmath.mp.prec = 100

#DESolver(g,f,D,XY=[],DXDY=[], basis="c", N=20)

plt.clf()
plt.grid()

z = [mpf('1.0'), mpf('1.6023737137301802297425736482307e-31'), mpf('-0.16517585839980632738648516797681'), mpf('-0.015814894942047426307588932087475'), mpf('0.058539305931879546102667602463209'), mpf('-0.029569829217583336132321942358925'), mpf('0.0085219452927753322032061615005972'), mpf('-0.0016866596617630518476504387265631'), mpf('0.00024687729746271588158312411225029'), mpf('-0.000027808702024518631237835582812343'), mpf('0.0000024703757462076574420308418689777'), mpf('-0.00000017586979899524848939896460884123'), mpf('0.000000010140592546484511633685960336243'), mpf('-4.7670738136802271113959246782061e-10'), mpf('1.8332442651735882451530036041177e-11'), mpf('-5.7696378030372220973333176725294e-13'), mpf('1.4822354505825177421674741460143e-14'), mpf('-3.0899657558382269427453793882812e-16'), mpf('5.1748832640204138945738231097496e-18'), mpf('-6.854972421120947725319456290931e-20'), mpf('7.0146238559603502604311569755064e-22'), mpf('-5.3449122501059695488023922472283e-24'), mpf('2.853579092145327252730116592938e-26'), mpf('-9.5204794276744854600980917554685e-29'), mpf('1.4932788959245979011758888993301e-31')]

D = [0,50]

n = 3
points = 200
dx = (D[1]-D[0])/(points-1)
xlist = [D[0]+i*dx for i in range(points)]


for i in range(15):

	# OUr ODE's Coefficients
	def g0(x):
		return n*polyeval(z,x)**(n-1)

	def g1(x):
		return 2/x

	def g2(x):
		return 1

	# What our ODE equals
	def f(x):

		return (n-1)*polyeval(z,x)**n

	g = [g0,g1,g2]
	yi1 = z.copy()
	z = DESolver(g,f,D=D,XY=[[0],[1]],DXDY=[[0],[0]],basis="c",N=25)
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

plt.show()
plt.clf()
plt.grid()
plt.plot(xlist,aylist,"r",label="$y=y_{25}$")
plt.title("Lane-Emden $n=3$ Domain $[0,50]$ using SDEA")
plt.legend()
print(z)
plt.show()
