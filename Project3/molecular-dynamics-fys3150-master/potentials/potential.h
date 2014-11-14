#pragma once
#include <string>
#include <vector>
#include <system.h>

class Potential
{
protected:
    double m_potentialEnergy;
    double m_pressure;
public:
    Potential();
    virtual ~Potential() {}
    virtual void calculateForces(System *system) = 0;
    double potentialEnergy();
    double pressure();
    void setPotentialEnergy(double potentialEnergy);
};
