#include <CellList.h>
#include <system.h>
#include <iostream>

using namespace std;


CellList::CellList()
{

}

CellList::~CellList()
{

}

void CellList::setup(System *system, double r_cut)
{
    nx = (system->systemSize().x()/r_cut);
    ny = (system->systemSize().y()/r_cut);
    nz = (system->systemSize().z()/r_cut);

    cout << "Creating cells: " << nx << ", " << ny << ", " << nz << endl;

    cells.resize(nx);
    for(int i = 0; i < nx; i++) {
        cells[i].resize(ny);
        for(int j = 0; j < nz; j++) {
            cells[i][j].resize(nz);
        }
    }
}

void CellList::update(System *system)
{
    clearList();

    for(int n = 0; n < system->atoms().size(); n++) {
        Atom* atom = system->atoms()[n];
        int i = (atom->position.x()/system->systemSize().x())*nx;
        int j = (atom->position.y()/system->systemSize().y())*ny;
        int k = (atom->position.z()/system->systemSize().z())*nz;
        //cells[i][j][k].push_back(atom);

        if(i< 0 || i >= nx || j< 0 || j >= ny || k< 0 || k >= nz) {
            cout << i << " , " << j << " , " << k << endl;
        }

        try
        {
            cells.at(i).at(j).at(k).push_back(atom);
        } catch(string exception) {
            cout << "Crashed because: " << exception << endl;
            cout << "Cell indices: " << i << " , " << j << " , " << k << endl;
            cout << "Atom position: " << atom->position << endl;
        }
    }
}

void CellList::clearList()
{
    for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
            for(int k = 0; k < nz; k++) {
                cells[i][j][k].clear();
            }
        }
    }
}
