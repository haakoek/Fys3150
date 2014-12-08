from numpy import *

def LennardJones(r,sigma):
	return 4.0*((sigma/r)**12 - (sigma/r)**6)
	
sigma = 3.405
N = 100
r = linspace(3,7,N+1)
u = zeros(N+1)

for i in xrange(0,N+1):
	u[i] = LennardJones(r[i],sigma)
	
from matplotlib import pyplot as plt

plt.plot(r,u)
plt.title('LennardJones Potential')
plt.xlabel('r_ij')
plt.ylabel('U(r_ij)/eps')
plt.savefig('LennardJones.png')
