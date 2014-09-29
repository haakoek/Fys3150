from numpy import *

def V(p):
	return p**2
	
p_min = 0
p_max = 5

p = linspace(p_min,p_max,100)
V_vec = zeros(len(p))

for i in range(1,len(p)-1):

	V_vec[i] = V(p[i])
	
V_vec[0] = max(V_vec)
V_vec[len(p)-1] = 1.01*max(V_vec)

from matplotlib import pyplot as plt

plt.plot(p[1:len(p)-2],V_vec[1:len(p)-2])
plt.suptitle('Harmonic oscillator potential with rho_max = %d' % (p_max), fontsize=20)
plt.xlabel('rho', fontsize=20)
plt.ylabel('V(rho)=rho**2', fontsize=20)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.savefig("harmonicpotential.png")

