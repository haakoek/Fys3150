from mdconfig import *
from numpy import *
from scitools.all import *

program = MD(dt=10, name="molecular-dynamics-fys3150") # dt in femtoseconds
numberOfUnitCells = [10]
temp = 500

for n in numberOfUnitCells:

	program.runNewSystem(numberOfUnitCells = n, FCCLatticeConstant = 5.26, temperature=temp)
	program.saveState(path="states/00_initialState")

	program.runThermalize(timesteps = 1001)
	program.saveState(path="states/02_thermalized")
	program.runThermostat(temperature = temp, timesteps = 1001)
	program.saveState(path="states/01_T%d" % temp)
	
	
	temperatureFile = open('tempratureFile.txt', 'r')
	temperature = []
	t = []
	
	for line in temperatureFile:
		tmp = line.split()
		t.append(float(tmp[0]))
		temperature.append(float(tmp[1]))
		
	temperatureFile.close()
		
	temperature = array(temperature)
	t = array(t)
	
	n_atoms = 4*n**3
	print n_atoms
	
	plot(t,temperature,show=False)
	legend('N(atoms)=%d' % n_atoms )
	xlabel('time[10^(-14) seconds]')
	ylabel('T[K]')
	hold('on')
	savefig('thermoStatPlot.png')
	
