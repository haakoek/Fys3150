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

using namespace std;

int main()
{
    double dt = UnitConverter::timeFromSI(1e-14); // You should try different values for dt as well.


    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    cout << "One unit of pressure is " << UnitConverter::pressureToSI(1.0) << " Pa" << endl;
    cout << "One unit of energy is " << UnitConverter::energyToSI(1.0) << " J" << endl;
    cout << "Sigma: " << UnitConverter::lengthFromAngstroms(3.405) << endl;


    double b = 5.26;
    double boltz = 1.308*10E-23;
    cout << 119.8*boltz << endl;
    int numberOfUnitCellsEachDimension = 5;
    //int numberOfAtoms = pow(numberOfUnitCellsEachDimension,3)*4;
    //vec3 systemSize(b*numberOfUnitCellsEachDimension,b*numberOfUnitCellsEachDimension, b*numberOfUnitCellsEachDimension);

    System system;
    double systemTemp = 1000.0;
    double sigma = UnitConverter::lengthFromAngstroms(3.405);
    system.m_rCut = 2.5*sigma;
    //system.setSystemSize(UnitConverter::lengthFromAngstroms(systemSize));
    system.createFCCLattice(numberOfUnitCellsEachDimension, UnitConverter::lengthFromAngstroms(b), systemTemp);

    system.setPotential(new LennardJones(sigma, 1.0)); // You must insert correct parameters here
    system.setIntegrator(new VelocityVerlet());



    //

    StatisticsSampler *statisticsSampler = new StatisticsSampler(); //
    statisticsSampler->sample(&system);
    statisticsSampler->sampleMomentum(&system);
    system.setSystemNetMomentum(statisticsSampler->netMomentum);
    system.removeMomentum();


    IO *movie = new IO(); // To write the state to file
    movie->open("movie.xyz");

    clock_t begin1 = clock();

    BerendsenThermostat myThermostat(1.2,500.0); // tau = 1.2, Tbath = 500.

    int numTimeSteps = 1000;
    for(int timestep=0; timestep<numTimeSteps; timestep++) {

        statisticsSampler->sample(&system);

        if(timestep % 100 == 0) {
            statisticsSampler->printSample(timestep);
        }

        //myThermostat.adjustVelocity(&system, dt, statisticsSampler->temperature);

        system.step(dt);

        movie->saveState(&system);
    }

    clock_t end1 = clock();

    double elapsed_secs1 = double(end1 - begin1) / CLOCKS_PER_SEC;

    cout << "*********************" << endl;
    cout << "Time Usage: " << elapsed_secs1 << endl;
    cout << "Number of atoms: " << system.atoms().size() << endl;
    cout << "kiloAtomTimesteps per second of atoms: " << (system.atoms().size()*numTimeSteps)/elapsed_secs1/1000 << endl;
    cout << "*********************" << endl;

    movie->close();

    return 0;
}
