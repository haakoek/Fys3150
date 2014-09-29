# Import libraries
from scitools.std import *
from os import system
"""
from os import system
from numpy import *
"""

# Compile c++ code first.
compile_string = "g++ main.cpp -o exec.out -O3 -llapack -lblas -larmadillo"
system(compile_string)

n_start = 40
n_end = 200 
rho_max = 15

n_step = linspace(n_start,n_end,((n_end - n_start)/10) + 1)

iterations = zeros(len(n_step))
computed_l1 = zeros(len(n_step))
computed_l2 = zeros(len(n_step))
computed_l3 = zeros(len(n_step))
time_vec = zeros(len(n_step))
time_vec_arma = zeros(len(n_step))

for i in xrange(len(n_step)):

	systemStr = "./exec.out version2 %d %d" % (n_step[i],rho_max)
	system(systemStr)
	inFile = open('version2_output.txt', 'r')
	line = inFile.read()
	tmp = line.split()
	iterations[i] = tmp[1]
	computed_l1[i] = tmp[2]
	computed_l2[i] = tmp[3]
	computed_l3[i] = tmp[4]
	time_vec[i] = tmp[5]
	time_vec_arma[i] = tmp[6]
	inFile.close()
	
from matplotlib import pyplot as plt
balle = 0
#Plot of number of similarity transforms before all non-diagonal elements are less than eps=1.0e-8
plt.figure(balle)
plt.plot(n_step,iterations)
plt.suptitle("n_step vs number of similarity transformations with rho_max=%d, eps=1.0e-8" % rho_max, fontsize=14)
plt.xlabel('n_step', fontsize=14)
plt.ylabel('number of similarity transforms', fontsize=14)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.savefig("nstep_iterations.png")
plt.hold('off')


plt.figure(balle+1)
plt.plot(n_step,time_vec)
plt.suptitle("Runtime of Jacobi's Rotation Algorithm, rho_max=%d eps=1.0e-8" % rho_max, fontsize=17)
plt.xlabel('n_step', fontsize=20)
plt.ylabel('time[s]', fontsize=20)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.savefig("Jacobitimeplot.png")


plt.figure(balle+2)
plt.plot(n_step,time_vec_arma)
plt.suptitle("Runtime of Armadillos eigensolver, rho_max=%d eps=1.0e-8" % rho_max, fontsize=17)
plt.xlabel('n_step', fontsize=20)
plt.ylabel('time[s]', fontsize=20)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.savefig("armatimeplot.png")

for i in xrange(len(n_step)):
	if n_step[i] >= 100:
		print "n_step=%d  l1c=%.5f  l2c=%.5f  l3c=%.5f" % (n_step[i],computed_l1[i], computed_l2[i], computed_l3[i])



