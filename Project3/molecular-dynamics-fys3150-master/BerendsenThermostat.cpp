#include <BerendsenThermostat.h>
#include <cmath>
#include <unitconverter.h>
#include <iostream>
#include <atom.h>
#include <math/vec3.h>

using namespace std;


BerendsenThermostat::BerendsenThermostat(double tau, double T_bath, double dt) :
    m_tau(tau),
    m_Tbath(T_bath),
    m_dt(dt)
{

}

BerendsenThermostat::~BerendsenThermostat()
{

}

void BerendsenThermostat::adjustVelocity(System* system, double systemTemp)
{
    double T = systemTemp;

    double gamma = sqrt(1.0 + (m_dt/m_tau)*((m_Tbath/T) - 1.0));


    for(int i = 0; i < system->atoms().size(); i++) {

        Atom *atom_i = system->atoms()[i];

        atom_i->velocity *= gamma;
    }
}
