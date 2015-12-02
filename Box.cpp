#include "Box.h"
double L;                 // size of whole cube
double a;               //LatticeA constant
double M;                //number of cubes in LatticeA
double rCutOff;     // cut-off on Lennard-Jones potential
double lambda =0;

Box::Box(int N, int steps, string flag) {
    Box::MCSteps = steps;
    Box::N = N;
    Box::Lattice = new double *[N];
    for (int i = 0; i < N; i++) {
        Box::Lattice[i] = new double[3];
    }

    Box::UAvg = new double [Box::MCSteps];
    Box::Var= new double [Box::MCSteps];

    // compute side of cube from number of particles and number density

    L = pow(N / rho, 1.0/3);
    rCutOff = 0.49 * L;
    // find M large enough to fit N atoms on an fcc LatticeA
    M = 1;
    while (4 * M * M * M < N)
        ++M;
    a = L / M;           // LatticeA constant of conventional cell
    double dx[4] = {0.0, 0.5, 0.5, 0.0};//
    double dy[4] = {0.0, 0.0, 0.5, 0.5};//
    double dz[4] = {0.0, 0.5, 0.0, 0.5};//

    int n = 0;                  // atoms placed so far
    for (int x = 0; x < M; x++){
        for (int y = 0; y < M; y++){
            for (int z = 0; z < M; z++){
                for (int p = 0; p < 4; p++){
                    if (n < N) {
                        if(flag.compare("Right") == 0) {
                            Lattice[n][0] = (0.25 + x + dx[p]) * a - 0.5 * L;
                            Lattice[n][1] = (0.25 + y + dy[p]) * a;
                            Lattice[n][2] = (0.25 + z + dz[p]) * a - 0.5 * L;
                        }else {

                            Lattice[n][0] = (0.25 + x + dx[p]) * a - 0.5 * L;
                            Lattice[n][1] = n == 0 ?
                                                 -(0.75 + y + dy[p]) * a :
                                                 -(0.25 + y + dy[p]) * a;
                            Lattice[n][2] = (0.25 + z + dz[p]) * a - 0.5 * L;
                        }
                        ++n;
                    }
                }
            }
        }
    }
}

double Box::MinImage(int i, int j) {
    // find separation using closest image convention
    double dr[3];
    double sum = 0;
    for (int d = 0; d < 3; d++) {
        dr[d] = Lattice[i][d] - Lattice[j][d];
        if (dr[d] >= 0.5*L) dr[d] -= L;
        if (dr[d] < -0.5*L) dr[d] += L;
        sum += dr[d] * dr[d];
    }
    return sum;
}

void Box::computeOrderParameter() {
    //compute the order parameter after each MonteCarlo time step
    lambda = 0;
    for (int i = 0; i < N; i++){
        lambda += cos((4 * PI * Lattice[i][0]) / a) +
                  cos((4 * PI * Lattice[i][1]) / a) +
                  cos((4 * PI * Lattice[i][2]) / a);
    }
    lambda /= -3 * (double) N;
}

void Box::MMSteps(int i, double d, string filename, double step) {

}

double Box::computeEnergy(string flag) {
    return 0;
}

void Box::writeHeaders(ofstream &f, int index, string *fnames) {

}



double Box::ran0(long *idum) {
    long     k;
    double   ans;
    *idum ^= MASK;
    k = (*idum)/IQ;
    *idum = IA*(*idum - k*IQ) - IR*k;
    if(*idum < 0) *idum += IM;
    ans=AM*(*idum);
    *idum ^= MASK;
    return ans;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK
#undef PI

