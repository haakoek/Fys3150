#include <potentials/lennardjones.h>
#include <cmath>
#include <iostream>

using namespace std;

LennardJones::LennardJones(double sigma, double epsilon) :
    m_sigma(sigma),
    m_epsilon(epsilon)
{

}

void LennardJones::calculateForces(System *system)
{
    for(int i = 0; i < system->atoms().size(); i++) {
        Atom *atom_i = system->atoms()[i];
        atom_i->force.setToZero();
    }

    m_potentialEnergy = 0; // Remember to compute this in the loop

    for(int i = 0; i < system->atoms().size(); i++) {
        for(int j = i+1; j < system->atoms().size(); j++) {

            Atom *atom_i = system->atoms()[i];
            Atom *atom_j = system->atoms()[j];

            vec3 r_ij = atom_i->position - atom_j->position;

            //Minimum Image Criterion

            if(fabs(r_ij.x()) > system->systemSize().x()*0.5) {
                if (r_ij.x() > 0) {
                    r_ij.setX(r_ij.x() - system->systemSize().x());
                } else {
                    r_ij.setX(r_ij.x() + system->systemSize().x());
                }
            }

            if(fabs(r_ij.y()) > system->systemSize().y()*0.5) {
                if (r_ij.y() > 0) {
                    r_ij.setY(r_ij.y() - system->systemSize().y());
                } else {
                    r_ij.setY(r_ij.y() + system->systemSize().y());
                }
            }

            if(fabs(r_ij.z()) > system->systemSize().z()*0.5) {
                if (r_ij.z() > 0) {
                    r_ij.setZ(r_ij.z() - system->systemSize().z());
                } else {
                    r_ij.setZ(r_ij.z() + system->systemSize().z());
                }
            }

            //End Minimum Image Criterion

            double dr = r_ij.length();


            double DU_dr = -(24.0*m_epsilon/dr)*(2.0*(pow(m_sigma/dr,12)) - (pow(m_sigma/dr,6)));

            double Fx = -DU_dr*(r_ij.x()/dr);
            double Fy = -DU_dr*(r_ij.y()/dr);
            double Fz = -DU_dr*(r_ij.z()/dr);

            vec3 Force_i(Fx,Fy,Fz);
            vec3 Force_j(-Fx,-Fy,-Fz);

            atom_i->force.add(Force_i);
            atom_j->force.add(Force_j);

            m_potentialEnergy = m_potentialEnergy + 4.0*m_epsilon*(pow(m_sigma,12)/pow(dr,12) - pow(m_sigma,6)/pow(dr,6));
        }
    }

    //cout << "PotEn:" << m_potentialEnergy << endl;
}
