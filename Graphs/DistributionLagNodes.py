from DESolver5 import *
from datetime import datetime
from main import *

mpmath.mp.dps = 300
mpmath.mp.prec = 300

plt.title(r"Distribution of the $30^{th}$ Laguerre nodes.")
plt.xlabel(r"$x$")

poly = LaguerreGen(31)
ps = rootfind(poly[30])
ps.sort()
for p in ps:
	plt.plot(p,0,"og",alpha=0.5)



plt.show()