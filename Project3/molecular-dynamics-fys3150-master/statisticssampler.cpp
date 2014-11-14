#include <statisticssampler.h>
#include <iostream>
#include <potentials/lennardjones.h>
#include <unitconverter.h>

using namespace std;

StatisticsSampler::StatisticsSampler(IO *myFileManager)
{
    temperature = 0.0;
    potentialEnergy = 0.0;
    kineticEnergy = 0.0;
    numberDensity = 0.0;
    pressure = 0.0;
    this->myFileManager = myFileManager;
}

StatisticsSampler::~StatisticsSampler()
{

}

void StatisticsSampler::sample(System *system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    sampleMomentum(system);
    sampleTemperature(system);
    sampleNumberDensity(system);
    potentialEnergy = system->potential()->potentialEnergy();
    pressure = numberDensity*temperature - (1.0/(3.0*system->volume()))*system->potential()->pressure();

    // write stuff to file

}

void StatisticsSampler::sampleKineticEnergy(System *system)
{
    double E_kinetic = 0;

    for(int i = 0; i < system->atoms().size();i++) {
        Atom *atom_i = system->atoms().at(i);
        double velocity_squared = atom_i->velocity.lengthSquared();
        E_kinetic = E_kinetic + 0.5*atom_i->mass()*velocity_squared;
    }

    kineticEnergy = E_kinetic;

    /*
    double boltz = 1.308*10E-23;
    double instantaneousTemperature = (2.0/3.0)*(E_kinetic/(system->atoms().size()*boltz));
    cout << "instTemp= " << instantaneousTemperature << endl;
    */
}

void StatisticsSampler::sampleTemperature(System *system) {
    temperature = (2.0/3.0)*(kineticEnergy/system->atoms().size());
}

vec3 StatisticsSampler::sampleMomentum(System *system) {
    netMomentum.setToZero();
    for(int i = 0; i < system->atoms().size(); i++) {
        Atom *atom = system->atoms().at(i);
        vec3 momentum = atom->velocity*atom->mass();
        netMomentum.add(momentum);
//        netMomentum.set(netMomentum.x() + atom->mass()*atom->velocity.x(),
//                        netMomentum.y() + atom->mass()*atom->velocity.y(),
//                        netMomentum.z() + atom->mass()*atom->velocity.z());
    }

    return netMomentum;
}

void StatisticsSampler::sampleNumberDensity(System *system) {
    numberDensity = (system->atoms().size())/(system->systemSize().x()*system->systemSize().y()*system->systemSize().z());
}

void StatisticsSampler::printSample(double timestep) {
    cout << "Timestep= " << timestep << endl;
    cout << "E_k= " << kineticEnergy << endl;
    cout << "E_p= " << potentialEnergy << endl;
    cout << "E_tot= "<< (kineticEnergy + potentialEnergy) << endl;
    cout << "Temperature= " << UnitConverter::temperatureToSI(temperature) << endl;
    cout << "**************" << endl;
}

void StatisticsSampler::writeStatisticsToFile(System* system, int timestep) {
    myFileManager->writeEnergyToFile("energyFile.txt", kineticEnergy, potentialEnergy, timestep);
    myFileManager->writetemperatureToFile("tempratureFile.txt", temperature, timestep);
    myFileManager->writePressureToFile(("pressureFile.txt"),pressure,timestep);
}

