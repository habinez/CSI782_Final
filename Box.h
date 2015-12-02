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
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

using namespace std;
extern long l = -1; // simulation parameters
extern const double rho = 1.0;         // density (number per unit volume)

class Box {
private:

public:
    static double **Lattice;     // positions,
    static int N; // number of particles
    static double  *UAvg;
    static double *Var;
    Box(int N, int steps, string flag);
    double ran0(long *);
    double MinImage(int i, int j);
    void computeOrderParameter();
    void MMSteps(int, double, string filename, double step);
    double computeEnergy(string flag);
    void writeHeaders(ofstream& f, int index, string *fnames);

    static int MCSteps;
};


#endif //FINAL_BOX_H
