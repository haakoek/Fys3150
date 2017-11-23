import numpy as np
import matplotlib.pyplot as plt
import sys

Nx = 10
Nt = 200
tfinal = 1.0

u = np.zeros(Nx+1)
u[Nx] = 1
x = np.linspace(0,1,Nx+1)

dt = tfinal/Nt
dx = 1.0/Nx
alpha = dt/dx**2

print "alpha: ", alpha
plt.figure(0)
plt.plot(x,u)
plt.savefig("U_t=0.000.png")
for j in range(0,Nt):
	t = (j+1)*dt
	u[1:Nx] = alpha*u[2:Nx+1] + (1.0-2.0*alpha)*u[1:Nx] + alpha*u[0:Nx-1]
	if( (j+1)%10 == 0 ):
		print t
		plt.figure(j+1)
		plt.plot(x,u)
		plt.savefig("U_t=%.3f.png" % t)



