import numpy as np
import matplotlib.pyplot as plt
import sys

Nx = 20
Nt = 1000
tfinal = 0.5

u = np.zeros(Nx+1)
u[Nx] = 1
x = np.linspace(0,1,Nx+1)

dt = tfinal/Nt
dx = 1.0/Nx
alpha = dt/dx**2

print "alpha: ", alpha
"""
plt.figure(0)
plt.plot(x,u)
plt.savefig("U_t=0.000.png")

for j in range(0,Nt):
	t = (j+1)*dt
	u[1:Nx] = alpha*u[2:Nx+1] + (1.0-2.0*alpha)*u[1:Nx] + alpha*u[0:Nx-1]
	if( (j+1)%50 == 0 ):
		print t
		plt.figure(j+1)
		plt.plot(x,u)
		plt.savefig("U_t=%.3f.png" % t)
"""
#2d equation

Nx = 40
Ny = 40
Nt = 4000
tfinal = 0.3

dt = tfinal/Nt
dx = 1.0/Nx
alpha = dt/dx**2

x    = np.linspace(0,1,Nx+1)
y    = np.linspace(0,1,Ny+1)
X, Y = np.meshgrid(x, y)

u = np.zeros((Nx+1,Ny+1))
u[:,0]    = 1
u[0,:]    = 0
u[Nx,:]   = 0
u[:,Ny]   = 2


plt.pcolor(X,Y,u)
plt.colorbar()
plt.savefig("U2d_t=0.000.png")

for j in range(0,Nt):
	t = (j+1)*dt
	u[1:Nx,1:Ny] = u[1:Nx,1:Ny] + alpha*(u[2:Nx+1,1:Ny] + u[1:Nx,2:Ny+1] - 4*u[1:Nx,1:Ny] + u[1:Nx,0:Ny-1] + u[0:Nx-1,1:Ny])
	if( (j+1)%50 == 0 ):
		print t
		plt.figure(j+1)
		plt.pcolor(X,Y,u)
		plt.colorbar()
		plt.savefig("U2d_t=%.3f.png" % t)
		plt.close()




