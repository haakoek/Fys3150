from numpy import *


energyFile = open('energyFile.txt', 'r')
tempratureFile = open('tempratureFile.txt', 'r')
pressureFile = open('pressureFile.txt','r')

kineticEnergy = []
potentialEnergy = []
temprature = []
pressure = []
time_vec = []

for line in energyFile:
	tmp = line.split()
	time_vec.append(float(tmp[0]))
	kineticEnergy.append(float(tmp[1]))
	potentialEnergy.append(float(tmp[2]))
	
for line in tempratureFile:
	tmp = line.split()
	temprature.append(float(tmp[1]))
	
for line in pressureFile:
	tmp = line.split()
	pressure.append(float(tmp[1]))

pressureFile.close()	
energyFile.close()
tempratureFile.close()
	
kineticEnergy = array(kineticEnergy)
potentialEnergy = array(potentialEnergy)
temprature = array(temprature)
pressure = array(pressure)
time_vec = array(time_vec)

print pressure

totalEnergy = zeros(len(kineticEnergy))

for i in range(0,len(kineticEnergy)):
	totalEnergy[i] = kineticEnergy[i] + potentialEnergy[i]
	
from scitools.all import *

figure()
plot(time_vec,kineticEnergy, '-r', time_vec, potentialEnergy, '-b', time_vec, totalEnergy, '-g',show=True)
xlabel('time[10^(-14) seconds]')
ylabel('Energy[units of energy]')
legend('kineticEnergy', 'potentialEnergy', 'totalEnergy')
savefig('energyPlot.png')

figure()
plot(time_vec, temprature,show=True)
xlabel('time[10^(-14) seconds]')
ylabel('T[K]')
savefig('temperaturePlot.png')


figure()
plot(time_vec,pressure,show=True)
xlabel('time[10^(-14) seconds]')
ylabel('P[units of pressure]')
savefig('pressurePlot.png')
