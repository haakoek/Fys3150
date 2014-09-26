# Import libraries
from scitools.std import *
from os import system
"""
from math import exp
import matplotlib.pyplot as plt
from os import system
from numpy import *
"""

# Compile c++ code first.
#system("g++ code.cpp -o code.out -O3")

nStart = 100
n      = 5000
sample = 10
N      = linspace(nStart, n, sample)
print N
rel_error  = zeros(len(N))
run_time_solver = zeros(len(N))
run_time_gauss = zeros(len(N))
run_time_LU = zeros(len(N))
run_time_random = zeros(len(N))

for i in xrange(len(N)) :

    # Call c++ code 
    systemStr = "./exec.out %d output1.txt" % N[i]
    system(systemStr)
    
    # Get result from file.
    inFile = open('output1.txt', 'r')
    line = inFile.read()
    tmp = line.split()
    rel_error[i] = tmp[0]
    run_time_solver[i] = tmp[1]
    print tmp
    if (n <= 10000):
    	run_time_gauss[i] = tmp[2]
    	run_time_LU[i] = tmp[3]
    	run_time_random[i] = tmp[4]
    
    inFile.close()
    

if(n>10000):
	figure()
	plot(N,rel_error,title='Maximum relative error vs. log(N), with n=%d sample points' % sample,
	 	xlabel='N',ylabel='eps(N)',
	 	log='x',
	 	axis=[nStart,n,0,-10],
	 	savefig="relerror.png")

	

if(n<=10000):
	hold('on')
	figure()
	plot(N,run_time_gauss, '-b', 
		 N, run_time_LU, '-g',
		 xlabel='N',
		 ylabel='t',
		 legend=('Gauss time','LU time'),
		 savefig="GaussLUTime.pdf",
		 show=False,
		 title="Runtime of armadillo LU and Gauss, with n=%d sample points" % sample)

	hold('on')

	figure()
	plot(N,run_time_solver,'-r',
	 	title='Runtime of solver, with n=%d sample points' % sample,
	 	xlabel='N',ylabel='t',
	 	savefig="rtsolver.png",
	 	show=False)
	
	hold('on') 	
	
	figure()
	plot(N,run_time_random, '-r',
		title='Runtime of Gauss with random matrix C, with n=%d sample points' % sample,
		xlabel='N', ylabel='t',
		savefig="random.png",
		show=False)
    
    
    
