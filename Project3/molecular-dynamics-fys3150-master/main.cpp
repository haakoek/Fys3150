#include <iostream>
#include <math/random.h>
#include <cmath>

#include <potentials/lennardjones.h>
#include <integrators/eulercromer.h>
#include <integrators/velocityverlet.h>
#include <system.h>
#include <statisticssampler.h>
#include <atom.h>
#include <io.h>
#include <unitconverter.h>
#include <time.h>
#include <BerendsenThermostat.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int main(int argc, char* argv[])
{
    double dt = UnitConverter::timeFromSI(1e-15); // You should try different values for dt as well.
    int numTimeSteps = 2000;
    int numberOfUnitCells = 5;
    double latticeConstant = 5.26;
    bool loadState = true;
    bool thermostatEnabled = false;
    double temperature = 200.0;
    double relaxationTime = dt*100;

    if(argc > 1) {
        dt = UnitConverter::timeFromSI(atof(argv[1])*1e-15);
        numTimeSteps = atoi(argv[2]);
        numberOfUnitCells = atoi(argv[3]);
        latticeConstant = atof(argv[4]);
        loadState = atoi(argv[5]);
        thermostatEnabled = atoi(argv[6]);
        temperature = atof(argv[7]);
    }

    System system;
    IO fileManager;


    double sigma = UnitConverter::lengthFromAngstroms(3.405);

    system.m_rCut = 2.5*sigma;

    if(loadState) {
        cout << "LOADED!!!" << endl;
        fileManager.load("state.bin",&system);
    } else {

        system.createFCCLattice(numberOfUnitCells, UnitConverter::lengthFromAngstroms(latticeConstant), temperature);
    }

    LennardJones* potential = new LennardJones(sigma,1.0);
    potential->useCellLists = true;
    system.setPotential(potential);
    system.setIntegrator(new VelocityVerlet());

    StatisticsSampler *statisticsSampler = new StatisticsSampler(&fileManager); //
    statisticsSampler->sample(&system);
    statisticsSampler->sampleMomentum(&system);
    system.setSystemNetMomentum(statisticsSampler->netMomentum);
    system.removeMomentum();

    IO *movie = new IO(); // To write the state to file
    movie->open("movie.xyz");

    clock_t begin1 = clock();

    BerendsenThermostat myThermostat(relaxationTime,UnitConverter::temperatureFromSI(temperature),dt);

    for(int timestep=0; timestep<numTimeSteps; timestep++) {

        system.step(dt);
        statisticsSampler->sample(&system);


        if(timestep % 200 == 0) {
            statisticsSampler->printSample(timestep);
            statisticsSampler->writeStatisticsToFile(&system,timestep);

        }


        if(thermostatEnabled) {
            myThermostat.adjustVelocity(&system, statisticsSampler->temperature);
        }

        movie->saveState(&system);

    }

    clock_t end1 = clock();

    double elapsed_secs1 = double(end1 - begin1) / CLOCKS_PER_SEC;
    double kiloAtomTimesteps = (system.atoms().size()*numTimeSteps)/elapsed_secs1/1000;
    int numberOfAtoms = system.atoms().size();

    /*
    cout << "*********************" << endl;
    cout << "Time Usage: " << elapsed_secs1 << endl;
    cout << "Number of atoms: " << system.atoms().size() << endl;
    cout << "kiloAtomTimesteps per second of atoms: " << (system.atoms().size()*numTimeSteps)/elapsed_secs1/1000 << endl;
    cout << "*********************" << endl;
    */

    movie->close();


    fileManager.save("state.bin",&system);

    return 0;
}
