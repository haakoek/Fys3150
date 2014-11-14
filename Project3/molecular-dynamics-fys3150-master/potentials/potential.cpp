#include <potentials/potential.h>

Potential::Potential() :
    m_potentialEnergy(0),
    m_pressure(0)
{

}

double Potential::potentialEnergy()
{
    return m_potentialEnergy;
}

double Potential::pressure()
{
    return m_pressure;
}

void Potential::setPotentialEnergy(double potentialEnergy)
{
    m_potentialEnergy = potentialEnergy;
}
