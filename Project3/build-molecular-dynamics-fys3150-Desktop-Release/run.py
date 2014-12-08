from mdconfig import *

program = MD(dt=10, name="molecular-dynamics-fys3150") # dt in femtoseconds
# Create new FCC lattice
numberOfUnitCells = 1
temperature = 500
program.runNewSystem(numberOfUnitCells = numberOfUnitCells, FCCLatticeConstant = 5.26, temperature=temperature) # lattice constant in angstroms
program.saveState(path="states/00_initialState")
#Run with thermostat for some time.
program.runThermostat(temperature = temperature, timesteps = 1001)
program.saveState(path="states/01_T%d" % temperature)
#Thermalize for some time (Note: thermalizing means without the thermostat, equilibrate)
program.runThermalize(timesteps = 2001)
program.saveState(path="states/02_thermalized")
#Sample
program.runThermalize(timesteps = 3001)

#Plotting
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

N = len(kineticEnergy)

totalEnergy = zeros(N)

for i in range(0,N):
	totalEnergy[i] = kineticEnergy[i] + potentialEnergy[i]
	
#Calculate heat capacity

expectedE = 0.0
expectedE2 = 0.0
Tavg = 0.0
numberOfAtoms = numberOfUnitCells**3*4
kb = 1.3806*10**(-23.0)
#kb = 1.0

for j in xrange(0,N):
	expectedE2 += kineticEnergy[j]**2
	expectedE  += kineticEnergy[j]
	Tavg += temprature[j]

expectedE2 = (1.0/N)*expectedE2
expectedE = (1.0/N)*expectedE
Tavg = (1.0/N)*Tavg

varEkinetic = expectedE2 - expectedE**2

Cv = (3.0*kb/2.0)*(1.0/(1.0-(2*N/(3.0*kb*Tavg**2))*(varEkinetic)))
m_arg = 6.63352088e-26
cv = Cv/m_arg

print "cv=%g" % cv
print "Cv=%g" % Cv
print "Pressure=%g" % (sum(pressure)/len(pressure))

"""
from scitools.all import *

figure()
plot(time_vec,kineticEnergy, '-r', time_vec, potentialEnergy, '-b', time_vec, totalEnergy, '-g',show=False)
title('Energyplot')
legend('kineticEnergy', 'potentialEnergy', 'totalEnergy')
savefig('energyPlot.png')

figure()
plot(time_vec, temprature,show=False)
title('Temperature Plot')
savefig('temperatureplot.png')

figure()
plot(time_vec,pressure,show=False)
title('Pressure Plot')
savefig('pressurePlot.png')
"""
#http://catalog.conveyorspneumatic.com/Asset/FLS%20Specific%20Heat%20Capacities%20of%20Gases.pdf
#program.saveState(path="states/03_statistics",plot=True)

# In your main function in your c++ program you can read in these variables as arguments like this: 
# First, we assume some default values, then if arguments are given, use them instead

# int numTimeSteps = 1000;
# double dt = UnitConverter::timeFromSI(1e-14); 
# int numUnitCells = 10;
# float latticeConstant = 5.26;
# bool loadState = false;
# bool thermostatEnabled = false;
# float temperature = 300;
# if(args>1) {
#     dt = UnitConverter::timeFromSI(atof(argv[1])*1e-15);
#     numTimeSteps = atoi(argv[2]);
#     numberOfUnitCells = atoi(argv[3]);
#     latticeConstant = atof(argv[4]);
#     loadState = atoi(argv[5]);
#     thermostatEnabled = atoi(argv[6]);
#     temperature = atof(argv[7]);
# }
