//#pragma once
#include <system.h>
#include <math/vec3.h>

class BerendsenThermostat
{
private:
    double m_tau;
    double m_Tbath;
    //int enabled = 1;
public:
    BerendsenThermostat(double tau, double T_bath);
    ~BerendsenThermostat();
    void adjustVelocity(System* system, double dt, double systemTemp);
};
