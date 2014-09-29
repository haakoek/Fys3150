# Import libraries
from matplotlib import pyplot as plt
from os import system
from numpy import *
"""
from os import system
from numpy import *
"""

# Compile c++ code first.
compile_string = "g++ main.cpp -o exec.out -O3 -llapack -lblas -larmadillo"
system(compile_string)

n = 200
rho_max = 5

systemStr = "./exec.out version2 %d %d" % (n,rho_max)
system(systemStr)

eigenvectors = open('eigenvectors.txt', 'r')

	
ground_state = zeros(n+1)
first_ex = zeros(n+1)
second_ex = zeros(n+1)
x = linspace(0,rho_max,n+1)
k = 0
balle = 0
	
for line in eigenvectors:
	
	tmp = line.split()
		
	ground_state[k] = tmp[0]
	first_ex[k] = tmp[1]
	second_ex[k] = tmp[2]
	
	k = k+1 
	

for i in xrange(len(ground_state)):
	ground_state[i] = ground_state[i]**2
	first_ex[i] = first_ex[i]**2
	second_ex[i] = second_ex[i]**2 
	
	
	
	


eigenvectors.close()	
	
plt.figure(balle)
plt.plot(x,ground_state)
plt.suptitle('Ground_state one electron', fontsize=20)
plt.xlabel('rho', fontsize=20)
plt.ylabel('u(rho)**2', fontsize=20)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.savefig("ground_state_one.png")
plt.hold('off')
	
	
plt.figure(balle+1)
plt.plot(x,first_ex)
plt.suptitle('First Excited one electron', fontsize=20)
plt.xlabel('rho', fontsize=20)
plt.ylabel('u(rho)**2', fontsize=20)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.savefig("first_ex_one.png")
plt.hold('off')
	
plt.figure(balle+2)
plt.plot(x,second_ex)
plt.suptitle('Second Excited state one electron', fontsize=20)
plt.xlabel('rho', fontsize=20)
plt.ylabel('u(rho)**2', fontsize=20)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.savefig("second_ex_one.png")
plt.hold('off')

