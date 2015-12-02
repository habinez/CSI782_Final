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
    double *Move(int, double);
    void UnMove(int index, double *delta);

public:
    static double **Lattice;     // positions,
    static string position;
    static int N; // number of particles
    static double  *UAvg;
    static double *Var;
    Box(int, int, string);
    ~Box();
    double MinImage(int, int);
    void computeOrderParameter();
    void MMSteps(int, double, string filename, double step);
    double computeEnergy();
    void writeHeaders(ofstream& f, int index, string *fnames);

    static int MCSteps;
};


#endif //FINAL_BOX_H
