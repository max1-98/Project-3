from DESolver5 import *
import numpy as np

mpmath.mp.dps = 100
mpmath.mp.prec = 100

#DESolver(g,f,D,XY=[],DXDY=[], basis="c", N=20)

plt.clf()
plt.grid()

z = [mpf('1.0'), mpf('3.3896367021215351013785211789495e-32'), mpf('-0.16517585839980632738648516797681'), mpf('-0.015814894942047426307588932087844'), mpf('0.058539305931879546102667602463948'), mpf('-0.029569829217583336132321942359369'), mpf('0.0085219452927753322032061615007698'), mpf('-0.0016866596617630518476504387266047'), mpf('0.00024687729746271588158312411225799'), mpf('-0.000027808702024518631237835582813379'), mpf('0.000002470375746207657442030841869083'), mpf('-0.00000017586979899524848939896460884988'), mpf('0.000000010140592546484511633685960336807'), mpf('-4.7670738136802271113959246784927e-10'), mpf('1.8332442651735882451530036042371e-11'), mpf('-5.7696378030372220973333176729455e-13'), mpf('1.4822354505825177421674741461265e-14'), mpf('-3.0899657558382269427453793885299e-16'), mpf('5.174883264020413894573823110193e-18'), mpf('-6.8549724211209477253194562915468e-20'), mpf('7.014623855960350260431156976168e-22'), mpf('-5.3449122501059695488023922477503e-24'), mpf('2.8535790921453272527301165932275e-26'), mpf('-9.5204794276744854600980917564642e-29'), mpf('1.4932788959245979011758888994904e-31')]
D = [0,50]

n = 3
points = 100
dx = (D[1]-D[0])/(points-1)
xlist = [D[0]+i*dx for i in range(points)]


for i in range(1):

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
		#plt.plot(xlist,aylist,":",label=lab)
		poep = 1


aylist1 = []
for ay in aylist:
	aylist1.append(float(ay))

aylist1 = np.array(aylist1,dtype=float)

plt.clf()
plt.grid()
plt.ylabel("y")
plt.xlabel("x")
plt.plot(xlist,aylist1,"-r",label="$y=y_{25}$")
plt.plot(xlist,np.array([1]*100,dtype=float),"-g", label="$y=1$")

plt.fill_between(xlist,aylist1,[1], color='r', alpha=0.4)

plt.title("Lane-Emden $n=3$ Domain $[0,50]$ vs $y=1$")
plt.legend()
plt.show()
