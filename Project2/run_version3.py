# Import libraries
#from scitools.std import *
from matplotlib import pyplot as plt
from os import system
from numpy import *


# Compile c++ code first.
compile_string = "g++ main.cpp -o exec.out -O3 -llapack -lblas -larmadillo"
system(compile_string)

omega_r = [0.01, 0.5, 1, 5]
rho_max = 3
n = 200
balle = 0
for om in omega_r:

	systemStr = "./exec.out version3 %d %d %f" % (n,rho_max,om)
	system(systemStr)
	
	inFile = open('version3_output.txt', 'r')
	
	ground_state = zeros(n+1)
	first_ex = zeros(n+1)
	second_ex = zeros(n+1)
	x = linspace(0,rho_max,n+1)
	k = 0
	for line in inFile:
		tmp = line.split()
		
		ground_state[k] = tmp[0]
		first_ex[k] = tmp[1]
		second_ex[k] = tmp[2]
		k = k+1 
	
	for i in xrange(len(ground_state)):
		ground_state[i] = ground_state[i]**2
		first_ex[i] = first_ex[i]**2
		second_ex[i] = second_ex[i]**2 
	
	
	plt.figure(balle)
	plt.plot(x,ground_state)
	plt.suptitle('Ground_state omega = %.2f rho_max = %.1f' % (om, rho_max), fontsize=20)
	plt.xlabel('rho', fontsize=20)
	plt.ylabel('u(rho)**2', fontsize=20)
	plt.rc('xtick', labelsize=15)
	plt.rc('ytick', labelsize=15)
	plt.savefig("ground_state%.3f.png" % om)
	plt.hold('off')
	
	
	plt.figure(balle+1)
	plt.plot(x,first_ex)
	plt.suptitle('First Excited omega = %.2f rho_max = %.1f' % (om, rho_max), fontsize=20)
	plt.xlabel('rho', fontsize=20)
	plt.ylabel('u(rho)**2', fontsize=20)
	plt.rc('xtick', labelsize=15)
	plt.rc('ytick', labelsize=15)
	plt.savefig("first_ex%.3f.png" % om)
	plt.hold('off')
	
	plt.figure(balle+2)
	plt.plot(x,second_ex)
	plt.suptitle('Second Excited omega = %.2f rho_max = %.1f' % (om, rho_max), fontsize=20)
	plt.xlabel('rho', fontsize=20)
	plt.ylabel('u(rho)**2', fontsize=20)
	plt.rc('xtick', labelsize=15)
	plt.rc('ytick', labelsize=15)
	plt.savefig("second_ex%.3f.png" % om)
	plt.hold('off')
	
	balle = balle + 3
		
	"""
	figure()
	plot(x,ground_state,savefig="ground_state%.3f.png" % om, show=False,title="ground_state w = %.3f rho_max = %.1f" % (om,rho_max),fontsize=30, xlabel="p", ylabel="u(p)")
	figure()
	plot(x,first_ex,savefig="first_ex%.3f.png" % om, show=False,title="first_ex w = %.3f rho_max = %.1f" % (om, rho_max), xlabel="p", ylabel="u(p)")
	figure()
	plot(x,second_ex,savefig="second_ex%.3f.png" % om, show=False,title="second_ex w = %.3f rho_max = %.1f" % (om, rho_max), xlabel="p", ylabel="u(p)")
"""	
