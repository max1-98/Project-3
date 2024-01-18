from DESolver5 import *

mpmath.mp.dps = 100
mpmath.mp.prec = 100

#  This procedure goes as follows:
#  We use Chebyshev to find a good solution on [0,1]
#  Then we sub x=1/u into our ODE
#  This transforms the [1,infinity] domain to [1,0] and we use Chebyshev to solve for g(u)=y(1/u) then using u=1/x we find y(x)=g(1/x)

# NOTE: For finding derivatives we can't use polydiff(g,x) for y(x) as y(x) is not a polynomial anymore...

plt.clf()
# Number of basis polynomials, iteration for first section, second sections and in total
N = 25
it1 = 2
it2 = 7
it = it1+it2

# Lane Emden Equation's n
n = 3


# First we use Chebyshev to create a good approximation on [0,1]

D = [0,1]
D1 = [1/25,1]

# Initial Approximation ODE 1
z1 = [1]

# Total Domain and forming xlist for plotting
D2 = [0,1/D1[0]]
points = 200
dx = (D2[1]-D2[0])/(points-1)
xlist = [D2[0]+i*dx for i in range(points)]

# Final Approximation

z1f = [mpf('0.99999999999999999999999999999921'), mpf('-1.3018332087515533682937137865996e-31'), mpf('-0.16666666666666666666657494256213'), mpf('-4.8521272380125980163931716297125e-20'), mpf('0.025000000000000005123414586587808'), mpf('-2.390503313687567174125995144807e-16'), mpf('-0.0037698412698349296215734013180841'), mpf('-1.0836514435923298967953321069921e-13'), mpf('0.00056859935460786316387036643016805'), mpf('-1.1137237087255282037982773934647e-11'), mpf('-0.000085763315069118233029452771292491'), mpf('-3.677508640227400121129509634735e-10'), mpf('0.000012937482421384866135513066246953'), mpf('-0.0000000046314848521305787433549709443833'), mpf('-0.0000019394151431815919676739146221118'), mpf('-0.00000002415334140726501425493433387313'), mpf('0.00000033425015749595918140523591629346'), mpf('-0.000000053054712287911120537563711555033'), mpf('0.000000011660574326245200770760477507555'), mpf('-0.000000046252776353620324225070104197859'), mpf('0.000000035532113158848638203286286679609'), mpf('-0.000000012714309265522876208827615057188'), mpf('0.0000000023378872346237623152687568849528'), mpf('-1.7503860344181271808908154111506e-10'), mpf('-1.7617860251196649876880883898289e-12')] 
z2f = [mpf('-1.8368353942665783447095014745823'), mpf('235.7032875521038060851635911018'), mpf('-12414.179048630092343708200278494'), mpf('363011.53417144573650892574369953'), mpf('-6907512.822717021717693918794707'), mpf('93124911.809237200243621052791977'), mpf('-935588344.82216224773727173431046'), mpf('7238709083.8553606797301070607033'), mpf('-44122051211.924468451620528963155'), mpf('215335634805.01901604017767888951'), mpf('-851419039421.97123606333696805605'), mpf('2750318196695.2891204232038096164'), mpf('-7299614547089.3924670451853918968'), mpf('15970776677049.277643806983904479'), mpf('-28834543782608.991283594637881016'), mpf('42905050375818.061523657971200629'), mpf('-52414125037579.319549893156744258'), mpf('52209130650416.607243237845142114'), mpf('-41949278091904.50219320360432862'), mpf('26753274454952.784105207977053353'), mpf('-13220133548835.726060227223102567'), mpf('4877303089692.0447376191857637578'), mpf('-1263711234107.2724799360689016039'), mpf('205037229531.32503044529951995514'), mpf('-15668665169.693706860061301272419')]
def j(x):
	if x < 1:
		return polyeval(z1f, x)
	else:
		return polyeval(z2f, 1/x)

eylist = [j(x) for x in xlist]

a = polyeval(z1f,1)
b = -polyeval(polydiff(z1f.copy()),1)
z2 = [0,a-b+1,b-2,1]


for i in range(it):
	if i < it1:
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
		z1 = DESolver(g,f,D=D,XY=[[0],[1]],DXDY=[[0],[0]],basis="c",N=N)
	else: 
		# OUr ODE's Coefficients
		def g0(x):
			return n*polyeval(z2,x)**(n-1)

		def g1(x):
			return 0

		def g2(x):
			return x**4

		# What our ODE equals
		def f(x):
			return (n-1)*polyeval(z2,x)**n

		g = [g0,g1,g2]
		z2 = DESolver(g,f,D=D1,XY=[[1],[a]],DXDY=[[1],[b]],basis="c",N=N)


	def h(x):
		if x < 1:
			return polyeval(z1, x)
		else:
			return polyeval(z2, 1/x)

	ylist = [h(x) for x in xlist]

	plt.title("Iteration: "+str(i+1))
	plt.plot(xlist,ylist,"--g",label="Iterate")
	plt.plot(xlist,eylist, "-r",label="Final Approximate Solution")
	plt.xlabel("x")
	plt.ylabel("y")
	ax = plt.gca()
	ax.set_ylim([-1.2, 1.8])
	plt.legend(loc="upper left")
	plt.show()
	plt.clf()

