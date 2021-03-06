#include <potentials/lennardjones.h>
#include <cmath>
#include <iostream>
#include <CellList.h>

using namespace std;

LennardJones::LennardJones(double sigma, double epsilon) :
    m_sigma(sigma),
    m_epsilon(epsilon),
    useCellLists(true)
{

}

void LennardJones::calculateForces(System *system)
{

    if(useCellLists) {

        calculateForcesWithCellList(system);
        return;
    }

    for(int i = 0; i < system->atoms().size(); i++) {
        Atom *atom_i = system->atoms()[i];
        atom_i->force.setToZero();
    }

    m_potentialEnergy = 0; // Remember to compute this in the loop
    vec3 systemSize = system->systemSize();

    for(int i = 0; i < system->atoms().size(); i++) {
        for(int j = i+1; j < system->atoms().size(); j++) {

            Atom *atom_i = system->atoms()[i];
            Atom *atom_j = system->atoms()[j];

            vec3 r_ij = atom_i->position - atom_j->position;

            //Minimum Image Criterion
            for(int a=0; a<3; a++) {
                if(r_ij[a] > 0.5*systemSize[a]) r_ij[a] -= systemSize[a];
                else if(r_ij[a] < -0.5*systemSize[a]) r_ij[a] += systemSize[a];
            }

            //End Minimum Image Criterion

            double dr = sqrt(r_ij.dot(r_ij)); //r_ij.length();


            double DU_dr = -(24.0*m_epsilon/dr)*(2.0*(pow(m_sigma/dr,12)) - (pow(m_sigma/dr,6)));

            double Fx = -DU_dr*(r_ij.x()/dr);
            double Fy = -DU_dr*(r_ij.y()/dr);
            double Fz = -DU_dr*(r_ij.z()/dr);

            vec3 Force_i(Fx,Fy,Fz);
            vec3 Force_j(-Fx,-Fy,-Fz);

            atom_i->force.add(Force_i);
            atom_j->force.add(Force_j);

            m_potentialEnergy = m_potentialEnergy + 4.0*m_epsilon*(pow(m_sigma,12)/pow(dr,12) - pow(m_sigma,6)/pow(dr,6)); //-U(r_c)
        }
    }

    //cout << "PotEn:" << m_potentialEnergy << endl;
}

void LennardJones::calculateForcesWithCellList(System *system)
{

    system->myCellist.update(system);

    for(int i = 0; i < system->atoms().size(); i++) {
        Atom *atom_i = system->atoms()[i];
        atom_i->force.setToZero();
    }



    m_potentialEnergy = 0;
    m_pressure = 0;
    double oneOverRcut = 1.0/system->m_rCut;
    double oneOverRcut2 = oneOverRcut*oneOverRcut;
    double oneOverRcut6 = oneOverRcut2*oneOverRcut2*oneOverRcut2;
    double sigma6 = m_sigma*m_sigma*m_sigma*m_sigma*m_sigma*m_sigma;
    double potentialEnergyAtRcut = (4*m_epsilon*sigma6*oneOverRcut6*(sigma6*oneOverRcut6 - 1));

    vec3 systemSize = system->systemSize();
    for(int cx = 0; cx < system->myCellist.nx; cx++) {
        for(int cy = 0; cy < system->myCellist.ny; cy++) {
            for(int cz = 0; cz < system->myCellist.nz; cz++) {
                vector<Atom *> &celle1 = system->myCellist.cells[cx][cy][cz];
                int celle1Size = celle1.size();
                for(int dx = -1; dx <= 1; dx++) {
                    for(int dy = -1; dy <= 1; dy++) {
                        for(int dz = -1; dz <=1; dz++) {
                            //Her må det skje noe....
                            vector<Atom *> &celle2 = system->myCellist.cells[(cx + dx + system->myCellist.nx) % system->myCellist.nx][(cy + dy + system->myCellist.ny) % system->myCellist.ny][(cz + dz + system->myCellist.nz) % system->myCellist.nz];
                            int celle2Size = celle2.size();

                            for(int i = 0; i < celle1Size; i++) {
                                Atom* atom_i = celle1[i];

                                for(int j = 0; j < celle2Size; j++) {
                                    Atom* atom_j = celle2[j];

                                    if(atom_i->index() <= atom_j->index()) {
                                        continue;
                                    }

                                    vec3 r_ij = atom_i->position;
                                    r_ij.addAndMultiply(atom_j->position, -1);

                                    for(int a=0; a<3; a++) {
                                        if(r_ij[a] > 0.5*systemSize[a]) r_ij[a] -= systemSize[a];
                                        else if(r_ij[a] < -0.5*systemSize[a]) r_ij[a] += systemSize[a];
                                    }

                                    double dr2 = r_ij.lengthSquared(); //sqrt(r_ij.dot(r_ij)); //r_ij.length();

                                    if(dr2 <= system->m_rCut*system->m_rCut) {


                                        double oneOverDr2 = 1.0/dr2;

                                        double oneOverDr6 = oneOverDr2*oneOverDr2*oneOverDr2;


                                        //double DU_dr = -(24.0*m_epsilon/dr)*(2.0*(pow(m_sigma/dr,12)) - (pow(m_sigma/dr,6)));
                                        double force = -24*m_epsilon*sigma6*oneOverDr6*(2*sigma6*oneOverDr6 - 1)*oneOverDr2;

                                        atom_i->force.addAndMultiply(r_ij, -force);
                                        atom_j->force.addAndMultiply(r_ij, force);

                                        m_pressure += force*sqrt(dr2)*dr2;
                                        m_potentialEnergy += (4*m_epsilon*sigma6*oneOverDr6*(sigma6*oneOverDr6 - 1) - potentialEnergyAtRcut);

                                        /*
                                        double Fx = -DU_dr*(r_ij.x()/dr);
                                        double Fy = -DU_dr*(r_ij.y()/dr);
                                        double Fz = -DU_dr*(r_ij.z()/dr);

                                        vec3 Force_i(Fx,Fy,Fz);
                                        vec3 Force_j(-Fx,-Fy,-Fz);

                                        atom_i->force.add(Force_i);
                                        atom_j->force.add(Force_j);

                                        m_potentialEnergy = m_potentialEnergy + 4.0*m_epsilon*(pow(m_sigma,12)/pow(dr,12) - pow(m_sigma,6)/pow(dr,6)) - 4.0*m_epsilon*(pow(m_sigma,12)/pow((system->m_rCut),12) - pow(m_sigma,6)/pow(system->m_rCut,6));
                                        m_pressure = m_pressure + (r_ij.dot(Force_i));
                                        */
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


}
