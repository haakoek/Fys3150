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
rho_max = 5

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
    

from scitools.all import *

#Plot of number of similarity transforms before all non-diagonal elements are less than eps=1.0e-8
figure()
plot(n_step,iterations,title="n_step vs number of similarity transformations with rho_max=%d, eps=1.0e-8" % rho_max,
	savefig="nstep_iterations.png", show=False,xlabel='n_step',ylabel='number of similarity transforms')

figure()
plot(n_step, time_vec,title="Runtime of Jacobi's Rotation Algorithm, rho_max=%d eps=1.0e-8" % rho_max,
	savefig="Jacobitimeplot.png",show=False,xlabel='n_step',ylabel='time[s]')

figure()
plot(n_step, time_vec_arma,title="Runtime of Armadillos eigensolver, rho_max=%d eps=1.0e-8" % rho_max,
	savefig="armatimeplot.png",show=False,xlabel='n_step',ylabel='time[s]')

for i in xrange(len(n_step)):
	if n_step[i] >= 100:
		print "n_step=%d  l1c=%.5f  l2c=%.5f  l3c=%.5f" % (n_step[i],computed_l1[i], computed_l2[i], computed_l3[i])



