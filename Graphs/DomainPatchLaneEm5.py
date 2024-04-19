from DESolver5 import *
from datetime import datetime

mpmath.mp.dps = 250
mpmath.mp.prec = 250

M=5
n = M


plt.title(r"4 Chebyshev sub-patches and 15 rational Chebyshev sub-patches (Lane-Emden $n=5$ on $[0,100]$)")

start=datetime.now()

# Total Domain [0,L1]
D = [0,100]

## Patches

# Chebyshev
c = 4

# Rational Chebyshev 
rc = 15

# Domain Chebyshev
D1=[0,10]

# Domain rational Chebyshev
D2 = [D1[1],D[1]]

## Equally wide sub patches for Chebsyshev
#pw = 10/c
#bl = [i*pw for i in range(c+1)]

## Chebyshev patches
bl = points(c-1,0,10,basis="c")
bl.reverse()
bl.append(D1[1])
bl.insert(0,0)

## Rational Chebyshev patches

# Domain in rational space
Dr = [1/(1+D2[1]),1/(1+D2[0])]
# In 1/(1+x) plane
brk = points(rc-1,Dr[0],Dr[1],basis="c")
brk.sort()
brk.append(Dr[1])
brk.insert(0,Dr[0])
bk = []

for x in brk:
	bk.append((1-x)/x)

bk.sort()

## List of functions on each patch 
zp = []
zp1 = []


### Chebyshev Sub Patches
# Initial guess Chebyshev patch 1
z = [1]
b1 = 1
b2 = 0 

for j in range(c):

	## Finding each iterate

	# Creating the domain
	Dc = [bl[j],bl[j+1]]
	for i in range(7):

		def g0(x):
			return M*polyeval(z,x)**(M-1)

		def g1(x):
			return 2/x

		def g2(x):
			return 1

		# What our ODE equals
		def f(x):
			return (M-1)*polyeval(z,x)**(M)

		g = [g0,g1,g2]
		yi1 = z.copy()
		
		z = DESolver(g,f,D=Dc,XY=[[Dc[0]],[b1]],DXDY=[[Dc[0]],[b2]],basis="c",N=25)

	# Creating initial conditions for the next patch
	b1 = polyeval(z,Dc[1])
	b2 = polyeval(polydiff(z),Dc[1])

	# Adds z to the list of z's 
	zp.append(z.copy())

	# Creates next initial guess
	z = [b1-b2*Dc[1],b2]


### Rational Chebyshev Subpatches
N1 = 25
it1 = 8

# Initial guess first patch
z1 = [0,(1+Dc[1])*(2*b1+(1+Dc[1])*b2),-(1+Dc[1])**2*((1+Dc[1])*b2+b1)]

# Initial Conditions 
b1 = b1
b2 = -(1+Dc[1])**2*b2

for j in range(rc):

	Dc = [brk[rc-(j+1)],brk[rc-j]]
	for i in range(it1):

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
		z2 = DESolver(g,f,D=Dc,XY=[[Dc[1]],[b1]],DXDY=[[Dc[1]],[b2]],basis="c",N=N1)
		z1 = z2.copy()

	# Adding converged solution to the list of solutions
	zp1.append(z2.copy())

	# Start Point/Initial Conditions on the next subpatch
	b1 = polyeval(z2,Dc[0])
	b2 = polyeval(polydiff(z2),Dc[0])
	z1 = [0,(2*b1-b2*Dc[0])/Dc[0],(b2*Dc[0]-b1)/Dc[0]**2]
print("Solution Found in: ", datetime.now()-start, " seconds.")

### Error Calculations
start=datetime.now()
sc = 0
ec = 0
for i in range(c):
	def err(x):
		return abs(x*polyeval(polydiffn(zp[i],2),x)+2*polyeval(polydiff(zp[i]),x)+x*polyeval(zp[i],x)**M)

	def er(x):
		return abs(polyeval(zp[i],x)-1/sqrt(1+x**2/3))

	ec += integrate(er,bl[i],bl[i+1])
	sc += integrate(err,bl[i],bl[i+1])


scr = 0
ecr = 0
for i in range(rc):
	def err(x):
		return abs((x-1)*x**3*polyeval(polydiffn(zp1[i],2),x)+2*x**3*polyeval(polydiff(zp1[i]),x)+(x-1)/x*polyeval(zp1[i],x)**M)*1/(x)**2

	def er(x):
		return abs(polyeval(zp1[i],1/(x+1))-1/sqrt(1+x**2/3))

	scr += integrate(err,brk[rc-i-1],brk[rc-i])
	ecr += integrate(er,bk[i],bk[i+1])

s = sc+scr
print("Exact Error Chebyshev Patch: ",ec)
print("Exact Error Rational Chebyshev Patch: ", ecr)
print("Residual Error Chebyshev Patch: ",sc)
print("Residual Error Rational Chebshev Patch: ",scr)
print("Residual Error on [0,b]: ",s)
print("LMRE on [0,b]: ",-math.log(1/D[1]*s,10))


print("Errors Found in: ", datetime.now()-start, " seconds.")
### 
start=datetime.now()

points = 1000
dx = (D[1]-D[0])/(points-1)
xlist = [D[0]+i*dx for i in range(points)]
plt.plot(xlist,[1/(sqrt(1+x**2/3)) for x in xlist],"-r",label=r"$y=\frac{1}{\sqrt{1+x^2/3}}$")

for i in range(c):
	points = 100
	dx = (bl[i+1]-bl[i])/(points-1)
	xlist = [bl[i]+j*dx for j in range(points)]
	ylist = [polyeval(zp[i],x) for x in xlist]
	plt.plot(xlist,ylist,"--g")

for i in range(rc):
	points = 200
	dx = (bk[i+1]-bk[i])/(points-1)
	xlist = [bk[i]+j*dx for j in range(points)]
	ylist = [polyeval(zp1[i],1/(1+x)) for x in xlist]

	if i == rc-1:
		plt.plot(xlist,ylist,"--g",label=r"$y_{\text{approx}}$")
	else:
		plt.plot(xlist,ylist,"--g")

plt.xlabel("x")
plt.ylabel("y")
plt.legend()
ax = plt.gca()
ax.set_ylim([0, 1.1])
print("Plots made in: ", datetime.now()-start, " seconds.")
plt.show()
