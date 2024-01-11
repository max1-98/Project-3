from main import *
import matplotlib.pyplot as plt
plt.clf()
plt.grid()
polys = LaguerreGen(10)
plt.ylim((-1,1))
plt.xlabel("x")
plt.ylabel("y")
a = 0
b = 23
points = 500

dx = (b-a)/points
xlist = [a+dx*i for i in range(points+1)]
c = ["b","g","r","c","m","y","k"]
for i in range(1,10,3):
	s = str(i)
	plt.plot(xlist, [polyeval(polys[i],x) for x in xlist],"-"+c[int((i-1)/2)],label=r"$\phi_"+s+"(x)$")
plt.plot(xlist,[0]*(points+1),"--r",label="x-axis")
plt.legend()
plt.show()