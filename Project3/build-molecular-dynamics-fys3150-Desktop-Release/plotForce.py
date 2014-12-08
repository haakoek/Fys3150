from numpy import*
from matplotlib import pyplot as plt


sigma = 3.405

def gradU(r,sigma):
	return (24.0/r)*(-2*(sigma/r)**12 + (sigma/r)**6)
	
r = linspace(2.0**(1.0/6.0)*sigma,7*sigma,1000)

gradU_vec = zeros(len(r))
eps = 10**(-5)
indeks_max = 0
r_max = 2.0**(1.0/6.0) 

for i in xrange(0,len(r)):
	gradU_vec[i] = gradU(r[i],sigma)
	
for i in xrange(0,len(r)-1):
	if (gradU_vec[i+1] >= gradU_vec[i]):
		indeks_max = i
		r_max = r[indeks_max]
		
counter = indeks_max
ratio = gradU_vec[counter]/gradU_vec[indeks_max]
r_cut = r_max
grad_max = gradU_vec[indeks_max] 
eps = grad_max/50.0

#If the force at r is 50 times smaller than grad_max it will be negligible... 

while(ratio >= eps):
	counter = counter+1
	ratio = gradU_vec[counter]/gradU_vec[indeks_max]
	r_cut = r[counter]

print "ratio=%g, r_c = %g, r_c/sigma=%g" %( ratio, r_cut, r_cut/sigma)
		
plt.plot(r,gradU_vec)
plt.xlabel('r_ij')
plt.ylabel('F(r_ij)/eps')
plt.savefig("gradplot.png")
