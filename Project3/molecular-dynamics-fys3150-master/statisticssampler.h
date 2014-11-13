#pragma once
#include <system.h>
#include <io.h>

class StatisticsSampler
{


public:
    vec3 netMomentum;
    double temperature;
    double kineticEnergy;
    double potentialEnergy;
    double numberDensity;
    double pressure;
    IO* myFileManager;

    StatisticsSampler(IO* myFileManager);
    ~StatisticsSampler();
    void sample(System *system);
    void sampleKineticEnergy(System *system);
    void printSample(double timestep);
    void sampleNumberDensity(System *system);
    void sampleTemperature(System *system);
    void writeStatisticsToFile(System *system, int timestep);
    vec3 sampleMomentum(System *system);
};
