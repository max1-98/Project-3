import math
import matplotlib.pyplot as plt
from main import *
import mpmath
import latex
from mpmath import *

mpmath.mp.dps = 100
mpmath.mp.prec = 100

plt.title(r"A plot of $y_2$ vs the exact solution $\sqrt{1+x^2}$")
plt.grid()
plt.xlabel(r"$x$")
plt.ylabel(r'$y$')



a = -1
b = 1
points = 200
dx = (b-a)/points
xlist = [a+dx*i for i in range(points)]

y0 = [1,0,0.4797158088320493,0,-0.0633498498333404]
def y0p(x):

	return polyeval(y0,x)

y0list = [y0p(x) for x in xlist]


def e(x):
	return math.sqrt(1+x**2)
eylist = [e(x) for x in xlist]

y1 = [1,0,0.3660]

def y1p(x):

	return polyeval(y1,x)

y1list = [y1p(x) for x in xlist]

plt.plot(xlist,eylist, "-",c="r",label=r"$y=\sqrt{1+x^2}$")
plt.plot(xlist,y1list, "--", c="y", label=r"$y=y_2$")
#plt.plot(xlist,y0list, "--",c="g",label=r"Quartic Solution")




plt.legend(loc="upper center")

plt.show()