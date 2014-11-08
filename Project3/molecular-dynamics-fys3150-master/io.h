#pragma once
#include <fstream>
#include <string.h>
#include <math/vec3.h>


class System;
using std::ofstream;
using std::string;
class IO
{
private:
    ofstream file;
public:
    IO();
    ~IO();

    void saveState(System *system);
    void open(char *filename);
    void close();
    void save(string filename, System* system);
    void load(string filename, System* system);

};
