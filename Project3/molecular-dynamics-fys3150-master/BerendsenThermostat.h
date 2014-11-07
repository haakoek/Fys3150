//#pragma once
#include <system.h>
#include <math/vec3.h>

class BerendsenThermostat
{
private:
    double m_tau;
    double m_Tbath;
    double m_dt;
    //int enabled = 1;
public:
    BerendsenThermostat(double tau, double T_bath ,double dt);
    ~BerendsenThermostat();
    void adjustVelocity(System* system,  double systemTemp);
};
