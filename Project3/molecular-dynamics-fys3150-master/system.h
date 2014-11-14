#pragma once
#include <vector>
#include <atom.h>
#include <math/vec3.h>
#include <iostream>
#include <CellList.h>
#include <string.h>


class Potential; class Integrator;
using std::vector;
using std::string;


class System
{
private:
    vec3 m_systemSize;
    vec3 m_systemNetMomentum;
    vector<Atom*> m_atoms;
    Potential *m_potential;
    Integrator *m_integrator;
    double m_currentTime;

    int m_steps;


public:
    System();
    ~System();
    void resetForcesOnAllAtoms();
    void createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double systemTemp);
    void applyPeriodicBoundaryConditions();
    void removeMomentum();
    void calculateForces();
    void step(double dt);
    double volume();
    double m_rCut;
    CellList myCellist;

    // Setters and getters
    vector<Atom *> &atoms() { return m_atoms; }
    vec3 systemSize() { return m_systemSize; }
    void setSystemSize(vec3 systemSize) { m_systemSize = systemSize; myCellist.setup(this,m_rCut); }
    Potential *potential() { return m_potential; }
    void setPotential(Potential *potential) { m_potential = potential; }
    double currentTime() { return m_currentTime; }
    void setCurrentTime(double currentTime) { m_currentTime = currentTime; }
    Integrator *integrator() { return m_integrator; }
    void setIntegrator(Integrator *integrator) { m_integrator = integrator; }
    int steps() { return m_steps; }
    void setSteps(int steps) { m_steps = steps; }
    void setSystemNetMomentum(vec3 systemNetMomentum) { m_systemNetMomentum = systemNetMomentum; }
};
