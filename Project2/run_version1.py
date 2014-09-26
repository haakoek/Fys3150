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

rho_max = linspace(1,13,13)
l1 = zeros(len(rho_max))
l2 = zeros(len(rho_max))
l3 = zeros(len(rho_max))


for i in xrange(len(rho_max)):

	systemStr = "./exec.out version1 %d" % (rho_max[i])
	system(systemStr)
	inFile = open('version1_output.txt', 'r')
	line = inFile.read()
	tmp = line.split()
	l1[i] = tmp[0]
	l2[i] = tmp[1]
	l3[i] = tmp[2]
	inFile.close()
	
for i in xrange(len(rho_max)):
	print "rho_max=%d  l1c=%.5f  l2c=%.5f  l3c=%.5f" % (rho_max[i],l1[i], l2[i], l3[i])

