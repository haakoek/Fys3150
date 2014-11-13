#include <io.h>
#include <system.h>
#include <atom.h>
#include <unitconverter.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>


using std::endl;
using std::cout;
using std::string;
using std::ios;
using std::ifstream;
using namespace std;


IO::IO()
{

}

IO::~IO() {
    close();
}

void IO::open(char *filename) {
    if(file.is_open()) {
        std::cout << "<IO.cpp> Error, tried to open file " << filename << ", but some file is already open." << endl;
        exit(1);
    }

    file.open(filename);
}

void IO::close() {
    if(file.is_open()) {
        file.close();
    }
}

// This saves the current state to a file following the xyz-standard (see http://en.wikipedia.org/wiki/XYZ_file_format )
void IO::saveState(System *system)
{
    file << system->atoms().size() << endl;
    file << "The is an optional comment line that can be empty." << endl;
    for(int n=0; n<system->atoms().size(); n++) {
        Atom *atom = system->atoms()[n];
        file << "Ar " << UnitConverter::lengthToAngstroms(atom->position.x()) << " " << UnitConverter::lengthToAngstroms(atom->position.y()) << " " << UnitConverter::lengthToAngstroms(atom->position.z()) << endl;
    }
}

void IO::save(string filename, System* system)
{
    ofstream stateFile(filename.c_str(), ios::out | ios::binary);

    if (!stateFile.is_open()) {
        cout << "Could not open file " << filename << ". Aborting!" << endl;
        exit(1);
    }

    int numberOfPhaseSpaceCoordinates = 6*system->atoms().size();
    double *phaseSpace = new double[numberOfPhaseSpaceCoordinates];
    vec3 systemSize = system->systemSize();

    int phaseSpaceCounter = 0;

    for(int i = 0; i < system->atoms().size(); i++) {
        Atom* atom = system->atoms()[i];
        phaseSpace[phaseSpaceCounter++] = atom->position[0];
        phaseSpace[phaseSpaceCounter++] = atom->position[1];
        phaseSpace[phaseSpaceCounter++] = atom->position[2];
        phaseSpace[phaseSpaceCounter++] = atom->velocity[0];
        phaseSpace[phaseSpaceCounter++] = atom->velocity[1];
        phaseSpace[phaseSpaceCounter++] = atom->velocity[2];
    }

    int numberOfAtoms = system->atoms().size();
    stateFile.write(reinterpret_cast<char*>(&numberOfAtoms),sizeof(int));
    stateFile.write(reinterpret_cast<char*>(phaseSpace), numberOfPhaseSpaceCoordinates*sizeof(double));
    stateFile.write(reinterpret_cast<char*>(&systemSize[0]),3*sizeof(double));
    stateFile.close();
    delete phaseSpace;
}

void IO::load(string filename, System* system)
{
    ifstream stateFile(filename.c_str(), ios::in | ios::binary);

    if (!stateFile.is_open()) {
        cout << "Could not open file " << filename << ". Aborting!" << endl;
        exit(1);
    }

    int numberOfAtoms;
    stateFile.read(reinterpret_cast<char*>(&numberOfAtoms), sizeof(int));
    double *phaseSpace = new double[6*numberOfAtoms];
    stateFile.read(reinterpret_cast<char*>(phaseSpace), 6*numberOfAtoms*sizeof(double));
    vec3 systemSize;
    stateFile.read(reinterpret_cast<char*>(&systemSize[0]), 3*sizeof(double));
    system->setSystemSize(systemSize);

    int phaseSpaceCounter = 0;
    //double mass = 39.948;
    system->atoms().clear();

    for(int i = 0; i < numberOfAtoms; i++) {
        Atom* atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        atom->position = vec3(phaseSpace[phaseSpaceCounter++], phaseSpace[phaseSpaceCounter++], phaseSpace[phaseSpaceCounter++]);
        atom->velocity = vec3(phaseSpace[phaseSpaceCounter++], phaseSpace[phaseSpaceCounter++], phaseSpace[phaseSpaceCounter++]);
        system->atoms().push_back(atom);
    }
    stateFile.close();
    delete phaseSpace;

}

void IO::writeEnergyToFile(string filename, double kineticEnergy, double potentialEnergy, int timestep)
{
    if(!energyFile.is_open()) {
        energyFile.open(filename.c_str(),ios::out);
    }

    energyFile << setiosflags(ios::showpoint | ios::uppercase);
    energyFile << setw(15) << setprecision(8) << timestep;
    energyFile << setw(15) << setprecision(8) << kineticEnergy;
    energyFile << setw(15) << setprecision(8) << potentialEnergy << endl;

}

void IO::writetemperatureToFile(string filename, double temperature, int timestep)
{
    if(!temperatureFile.is_open()) {
        temperatureFile.open(filename.c_str(),ios::out);
    }

    temperatureFile << setiosflags(ios::showpoint | ios::uppercase);
    temperatureFile << setw(15) << setprecision(8) << timestep;
    temperatureFile << setw(15) << setprecision(8) << UnitConverter::temperatureToSI(temperature) << endl;
}
