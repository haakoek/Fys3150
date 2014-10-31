#pragma once
#include <system.h>

class StatisticsSampler
{


public:
    vec3 netMomentum;
    double temperature;
    double kineticEnergy;
    double potentialEnergy;

    StatisticsSampler();
    ~StatisticsSampler();
    void sample(System *system);
    void sampleKineticEnergy(System *system);
    void printSample(double timestep);
    void sampleTemperature(System *system);
    vec3 sampleMomentum(System *system);
};
