#include <integrators/velocityverlet.h>
#include <system.h>

VelocityVerlet::VelocityVerlet() :
    m_firstStep(true) // This will set the variable m_firstStep to false when the object is created
{

}

VelocityVerlet::~VelocityVerlet()
{

}

void VelocityVerlet::firstKick(System *system, double dt)
{
    m_firstStep = false;
    system->calculateForces();
    halfKick(system, dt);
}

void VelocityVerlet::halfKick(System *system, double dt)
{
    for(int n = 0; n < system->atoms().size(); n++) {
        Atom *atom = system->atoms()[n];
        double timeStepDividedByMass = dt/(atom->mass()*2);
        atom->velocity.addAndMultiply(atom->force, timeStepDividedByMass); // v(t+0.5dt) = v(t) + F/m*dt
    }
}

void VelocityVerlet::move(System *system, double dt)
{
    for(int n = 0; n < system->atoms().size(); n++) {
        Atom *atom = system->atoms()[n];
        atom->position.addAndMultiply(atom->velocity,dt); // r(t+dt) = r(t) + v(t+0.5dt)*dt
    }
}

void VelocityVerlet::integrate(System *system, double dt)
{
    if(m_firstStep) {
        firstKick(system, dt);
    } else {
        halfKick(system, dt);
    }

    move(system, dt);
    system->applyPeriodicBoundaryConditions();
    system->calculateForces();
    halfKick(system, dt);
}
