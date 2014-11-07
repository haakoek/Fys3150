#include <system.h>
#include <atom.h>
#include <integrators/integrator.h>
#include <potentials/potential.h>
#include <iostream>
#include <unitconverter.h>

using namespace std;

System::System() :
    m_potential(0),
    m_integrator(0),
    m_currentTime(0),
    m_steps(0),
    m_rCut(2.5*UnitConverter::lengthFromAngstroms(3.405))
{

}

System::~System()
{
    delete m_potential;
    delete m_integrator;
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
    for(int i = 0; i < m_atoms.size(); i++) {
        Atom *atom = m_atoms[i];

       if (atom->position.x() < 0) {
            atom->position.addX(m_systemSize.x());
       } else if(atom->position.x() >= m_systemSize.x()) {
            atom->position.addX(-m_systemSize.x());
       }

       if (atom->position.y() < 0) {
           atom->position.addY(m_systemSize.y());
       } else if(atom->position.y() >= m_systemSize.y()) {
           atom->position.addY(-m_systemSize.y());
       }

       if (atom->position.z() < 0) {
           atom->position.addZ(m_systemSize.z());
       } else if(atom->position.z() >= m_systemSize.z()) {
           atom->position.addZ(-m_systemSize.z());
       }
    }
}

void System::removeMomentum() {
    // Initially, when the atoms are given random velocities, there is a non-zero net momentum. We don't want any drift in the system, so we need to remove it.
    vec3 netMomentumPerAtom = m_systemNetMomentum/m_atoms.size();

    for(int i = 0; i < m_atoms.size(); i++) {
        Atom *atom = m_atoms[i];
        vec3 deltaVelocity = netMomentumPerAtom/atom->mass()*(-1);
        atom->velocity.add(deltaVelocity);
//        atom->velocity.setX(atom->velocity.x() - (m_systemNetMomentum.x()/(atom->mass()*m_atoms.size())));
//        atom->velocity.setY(atom->velocity.y() - (m_systemNetMomentum.y()/(atom->mass()*m_atoms.size())));
//        atom->velocity.setZ(atom->velocity.z() - (m_systemNetMomentum.z()/(atom->mass()*m_atoms.size())));
    }

}

void System::resetForcesOnAllAtoms() {

}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double systemTemp) {

    int N = numberOfUnitCellsEachDimension;
    double b = latticeConstant;

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            for(int k = 0; k < N; k++) {

                Atom * local_atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom * local_atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom * local_atom3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom * local_atom4 = new Atom(UnitConverter::massFromSI(6.63352088e-26));

                local_atom1->resetVelocityMaxwellian(UnitConverter::temperatureFromSI(systemTemp));
                local_atom2->resetVelocityMaxwellian(UnitConverter::temperatureFromSI(systemTemp));
                local_atom3->resetVelocityMaxwellian(UnitConverter::temperatureFromSI(systemTemp));
                local_atom4->resetVelocityMaxwellian(UnitConverter::temperatureFromSI(systemTemp));

                atoms().push_back(local_atom1);
                atoms().push_back(local_atom2);
                atoms().push_back(local_atom3);
                atoms().push_back(local_atom4);

                local_atom1->position.set(i*b,j*b,k*b);
                local_atom2->position.set(b*(0.5+i),b*(0.5+j),b*k);
                local_atom3->position.set(i*b,b*(0.5+j),b*(0.5+k));
                local_atom4->position.set(b*(0.5+i),j*b,(k+0.5)*b);
            }
        }

    }

    setSystemSize(vec3 (N*b,N*b,N*b));

}

void System::calculateForces() {
    resetForcesOnAllAtoms();
    m_potential->calculateForces(this);
}

void System::step(double dt) {
    //cout << "System sin step" << endl;
    m_integrator->integrate(this, dt);
    m_steps++;
    m_currentTime += dt;
}
