#include <statisticssampler.h>
#include <iostream>
#include <potentials/lennardjones.h>
#include <unitconverter.h>

using namespace std;

StatisticsSampler::StatisticsSampler()
{
    temperature = 0.0;
    potentialEnergy = 0.0;
    kineticEnergy = 0.0;
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
    potentialEnergy = system->potential()->potentialEnergy();
    // ...
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

void StatisticsSampler::printSample(double timestep) {
    cout << "Timestep= " << timestep << endl;
    cout << "E_k= " << kineticEnergy << endl;
    cout << "E_p= " << potentialEnergy << endl;
    cout << "E_tot= "<< (kineticEnergy + potentialEnergy) << endl;
    cout << "Temperature= " << UnitConverter::temperatureToSI(temperature) << endl;
    cout << "**************" << endl;
}
