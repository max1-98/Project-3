from mpmath import *
import matplotlib.pyplot as plt
import math 
from main import *
from DESolver5 import *
mpmath.mp.dps = 100
mpmath.mp.prec = 100

A = -0.9
B = 0.9
N = 3
p = polygen(N,A,B)
x = points(N-1,A,B)
print(x)
print(p)


iterations = 100
a = [0.5,1,0]



# Creates each row of f

# Checked
def f(a,i,p,x):

	if i == 0:
		m = 0

		for k in range(len(a)):
			m += a[k]*polyeval(p[k],0)

		return m

	s = 0
	t = 0
	# Evaluates yN at xi	
	for k in range(len(a)):
		s += polyeval(p[k],x[i-1])*a[k]

	# Evaluates yN' at xi
			
	for l in range(len(a)):
		t += polyeval(polydiff(p[l].copy()),x[i-1])*a[l]

	return cos(s)*t-1

# Creates the Vector f
def FU(a,p,x):
	F = []

	for i in range(len(a)):
		F.append(f(a,i,p,x))

	return F

# Creates each cell in the Jacobian
def J(a,i,p,x,s,t,j):
	# a is a list
	# p is a polynomial
	# x is a collocation point
	# i represents which a_i we are talking about

	# First equation is just a0phi_0(0)+a1phi_1(0)...
	if i == 0:
		return polyeval(p[j],0)


	return -polyeval(p[j],x)*sin(s)*t+cos(s)*polyeval(polydiff(p[j].copy()),x)

# Creates the Full Jacobian
def Jacobian(a,p,x):
	Ja = []
	for j in range(len(a)):
		R = []

		s = 0
		t = 0
		

		if j != 0:
			# Evaluates yN at xi	
			for k in range(len(a)):
				s += polyeval(p[k],x[j-1])*a[k]

			# Evaluates yN' at xi
			
			for l in range(len(a)):
				t += polyeval(polydiff(p[l].copy()),x[j-1])*a[l]

		for i in range(len(a)):
			if j == 0:
				R.append(J(a,j,p,0,s,t,i))
			else:
				R.append(J(a,j,p,x[j-1],s,t,i))
		Ja.append(R)

	return Ja


for i in range(iterations):

	
	JA = matrix(Jacobian(a,p,x))
	F = matrix(FU(a,p,x))
	a = matrix(a)

	c = JA*a-F
	a = mpmath.lu_solve(JA,c)



n = 200
dx = (B-A)/n
xlist = [A+dx*i for i in range(n)]

eylist = [math.asin(x) for x in xlist]

def af(x,p,a):

	m = 0
	for i in range(len(a)):
		m += polyeval(p[i],x)*a[i]

	return m
aylist = [af(x,p,a) for x in xlist]


plt.plot(xlist,eylist,"-r")
plt.plot(xlist,aylist,"-g")

def diffeq(x,p,a):
	s = 0
	for i in range(len(a)):
		s += polyeval(polydiff(p[i].copy()),x)*a[i]
	return cos(af(x,p,a))*s


dylist = [diffeq(x,p,a) for x in xlist]
plt.plot(xlist,dylist,"-y")


for xi in x:
	print(diffeq(xi,p,a))

print(a)

plt.show()








