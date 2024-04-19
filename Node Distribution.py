from DESolver5 import *

# Node Distribution
ps = points(3,-1,1,basis="c")
ps.sort()
ps.append(1)
ps.insert(0,-1)
print(ps)
pc = points(4*20,-1,1,basis="c")

w = 20
for i in range(len(pc)):
	
	if i ==0:
		plt.plot(pc[i],0.85,"ob",alpha=0.5,label="Nodes: No patching")
	else:
		plt.plot(pc[i],0.85,"ob",alpha=0.5)

for i in range(4):

	pl = points(w,ps[i],ps[i+1],basis="c")
	print(pl)
	pl1 = points(w,-1+i*1/2,-1+(i+1)*1/2,basis="c")
	for j in range(len(pl)):
		if i == 3 and j == 3:
			plt.plot(pl[j],1,"or",alpha=0.5,label="Nodes: Chebyshev Patches")
			plt.plot(pl1[j],0.5,"og",alpha=0.5,label="Nodes: Equal Width Patches")
		else:
			plt.plot(pl[j],1,"or", alpha=0.5)
			plt.plot(pl1[j],0.5,"og",alpha=0.5)

plt.title("Plot of the distribution of the Chebyshev nodes along x: Chebyshev patches vs equal width patches vs no patching")
plt.xlabel("x-value of a node")
plt.legend()
plt.show()

