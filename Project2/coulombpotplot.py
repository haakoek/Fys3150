from numpy import *
import sys

def V(p,omega):
	return omega**2*p**2 + 1.0/p
	
p_min = 0
p_max = int(sys.argv[1])
omega = float(sys.argv[2])

p = linspace(p_min,p_max,100)
V_vec = zeros(len(p))

for i in range(1,len(p)-1):

	V_vec[i] = V(p[i], omega)
	
V_vec[0] = max(V_vec)
V_vec[len(p)-1] = 1.01*max(V_vec)

from matplotlib import pyplot as plt

plt.plot(p[1:len(p)-2],V_vec[1:len(p)-2])
plt.suptitle('Coulomb potential with omega_r = %.2f, rho_max = %d' % (omega, p_max), fontsize=20)
plt.xlabel('rho', fontsize=20)
plt.ylabel('V(rho)=(omega_r**2)*rho**2 + 1/rho', fontsize=20)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.savefig("coulombpotential%.3f.png" % omega)


