import math
import matplotlib.pyplot as plt


def af2(x):
	return 1+0.47971580883204934*x**2 -0.06334984983334041*x**4

def af(x):

	return 1+0.36602540378443865*x**2

def ef(x):

	return math.sqrt(x**2+1)


points = 100
dx = 2/points
xlist = [-1+i*dx for i in range(points)]

aylist = [af(x) for x in xlist]
eylist = [ef(x) for x in xlist]
ay2list = [af2(x) for x in xlist]

plt.plot(xlist, aylist, "-g")
plt.plot(xlist, eylist, "-r")
plt.plot(xlist, ay2list, "-b")

plt.show()