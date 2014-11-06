#ifndef CELLLIST_H
#define CELLLIST_H



#include <atom.h>
#include <vector>

using std::vector;
class System;
class CellList
{
private:

public:
    CellList();
    ~CellList();
    vector<vector<vector<vector<Atom *> > > > cells;
    int nx;
    int ny;
    int nz;
    void setup(System *system, double r_cut);
    void update(System* system);
    void clearList();
};

#endif // CELLLIST_H
