#pragma once
#include <potentials/potential.h>



class LennardJones : public Potential
{
private:
    double m_sigma;
    double m_epsilon;

public:
    bool useCellLists;
    LennardJones(double sigma, double epsilon);
    ~LennardJones() {}
    virtual void calculateForces(System *system);
    void calculateForcesWithCellList(System *system);
};
