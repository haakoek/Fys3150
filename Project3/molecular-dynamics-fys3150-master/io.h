#pragma once
#include <fstream>
#include <string.h>
#include <math/vec3.h>


class System;
using std::ofstream;
using std::string;
class IO
{
private:
    ofstream file;

public:
    IO();
    ~IO();
    ofstream energyFile;
    ofstream temperatureFile;
    ofstream pressureFile;

    void saveState(System *system);
    void open(char *filename);
    void close();
    void save(string filename, System* system);
    void load(string filename, System* system);
    void writeEnergyToFile(string filename, double kineticEnergy, double potentialEnergy, int timestep);
    void writetemperatureToFile(string filename, double temperature, int timestep);
    void writePressureToFile(string filename, System* system);

};
