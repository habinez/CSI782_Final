#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>

#ifndef FINAL_BOX_H
#define FINAL_BOX_H
# define PI           3.14159265358979323846  /* pi */

using namespace std;


class Box {
private:
    vector<double> Move(int, double);
    void UnMove(int index, vector<double> delta);

public:
    double **Lattice;     // positions,
    string position;
    int N; // number of particles
    double  *UAvg;
    double *Var;
    Box(int, int, string);
    ~Box();
    double MinImage(int, int);
    void computeOrderParameter();
    void MMC(int, double, string filename, double step);
    void MMC(int, double, string filename, double step, Box B);
    double computeEnergy();
    void writeHeaders(ofstream& f, int index, string *fnames);

    int MCSteps;
};


#endif //FINAL_BOX_H
