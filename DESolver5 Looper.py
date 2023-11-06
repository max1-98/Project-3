from DESolver5 import *

mpmath.mp.dps = 500
mpmath.mp.prec = 500

#DESolver(g,f,D,XY=[],DXDY=[], basis="c", N=20)

def g0(x):
	return x*2*polyeval(z,x)

def g1(x):
	return 2

def g2(x):
	return x

# What our ODE equals
def f(x):

	return x*1*polyeval(z,x)**2



z = [1,0,-1/6,0,1/60]
D = [0,10]
g = [g0,g1,g2]

points = 100
dx = (D[1]-D[0])/(points-1)
xlist = [D[0]+i*dx for i in range(points)]
eylist = [polyeval([mpf('1.0000000000000000000000000000047'), mpf('-2.2602706304164754708824906134263e-31'), mpf('-0.16666518468218911239248943440978'), mpf('-0.000047835981503579706215130103926251'), mpf('0.01696966729136715209864594789421'), mpf('-0.00082768565738339960510628879823089'), mpf('-0.00021315033106398116541657630397451'), mpf('-0.001142211121455447520625057207395'), mpf('0.00079340538924191988255744009223085'), mpf('-0.00025719796029884002148253633284648'), mpf('0.000049292913255131130572971690215288'), mpf('-0.0000055123584567614088154532744596407'), mpf('0.00000021699557462006678757781517297784'), mpf('0.000000035770139645078675023748106876978'), mpf('-0.0000000071864691730238485756714843626306'), mpf('6.3167756019371508807520732423845e-10'), mpf('-3.1661312697384011647279134005873e-11'), mpf('8.536252016481749644899324395377e-13'), mpf('-8.0781727545594862508031761285327e-15'), mpf('-6.7176119461117863016612982273469e-17')],x) for x in xlist]



for i in range(1):

	def g1(x):
		return 2

	def g2(x):
		return x

	# What our ODE equals
	def f(x):

		return x*1*polyeval(z,x)**2

	z = DESolver(g,f,D=D,XY=[[0],[1]],DXDY=[[0],[0]],basis="c",N=30)
	print(i)
	print("This is solution #", i, ": ",z)
	aylist = [polyeval(z,x) for x in xlist]

	plt.plot(xlist,eylist,"r",label="20 OP")
	plt.plot(xlist,aylist,"--g",label="50 OP")
	plt.show()
	plt.clf()


print(z)





