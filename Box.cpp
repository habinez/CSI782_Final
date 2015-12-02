#include "Box.h"
double L;                 // size of whole cube
double a;               //LatticeA constant
double M;                //number of cubes in LatticeA
double rCutOff;     // cut-off on Lennard-Jones potential
double lambda =0;
double ran0(long *idum);
long l = -1; // simulation parameters
const double rho = 1.0;         // density (number per unit volume)

Box::Box(int N, int steps, string flag) {
    Box::MCSteps = steps;
    Box::N = N;
    Box::position = flag;
    Box::Lattice = new double *[N];
    for (int i = 0; i < N; i++) {
        Box::Lattice[i] = new double[3];
    }

    Box::UAvg = new double [Box::MCSteps];
    Box::Var= new double [Box::MCSteps];

    // compute side of cube from number of particles and number density

    L = pow(N / rho, 1.0/3);
    rCutOff = 0.49 * L;
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
                        if(position.compare("Right") == 0) {
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

void Box::MMSteps(int s, double T, string filename, double step) {
    ofstream file;
    file.precision(8);
    file.open(filename, ios::app);

    long ll = -1;
    double Ui = 0;
    double U = 0;
    double delta_r = step; //0.03;
    double AcceptanceRatio = 0;
    computeOrderParameter();

    vector<int> indexes;
    for (int i = 1; i < N; ++i) indexes.push_back(i);

    while (indexes.size() != 0) {
        // choose a random point, here we are using idx as index
        random_shuffle(indexes.begin(), indexes.end());
        int idx = indexes[indexes.size() - 1];
        indexes.pop_back();
        double Uold = computeEnergy();

        double *delta= Move(idx, delta_r);
        // Apply periodic boundary conditions
        for (int k = 0; k < 3; k++){
            if (Lattice[idx][k] < 0) Lattice[idx][k] += L;
            if (Lattice[idx][k] >= L)Lattice[idx][k] -= L;
            if (Lattice[idx][k] < 0) Lattice[idx][k] += L;
            if (Lattice[idx][k] >= L)Lattice[idx][k] -= L;

        }

        //check energy difference due to move
        double Unew = computeEnergy();
        double delta_U = Unew - Uold;
        if (delta_U > 0 && ran0(&ll) > exp(-(delta_U / T))) { //move is rejected
            UnMove(idx, delta);
            Ui = Uold;
        }
        else {//move is accepted
            Ui = Unew;
            AcceptanceRatio +=1;
        }
        U += (Ui)/((double)N);
        delete delta;
    }

    //compute the average energy and standard deviation using the method on
    UAvg[s] = s == 0 ? U : UAvg[s - 1] + ((U - UAvg[s - 1]) / s);

    Var[s] = s == 0 ? 0 : Var[s - 1] * s + (U - UAvg[s - 1]) * (Ui - UAvg[s]);
    Var[s] /= (s + 1);
    computeOrderParameter();
    AcceptanceRatio = (AcceptanceRatio * 100)/(double)N;

    file << s << "\t " << setw(12) << U << "\t            ";
    file << setw(12) << UAvg[s] << "\t  ";//running average
    file << setw(12) << Var[s] << "\t  ";//running variance
    file << setw(12) << lambda << "\t ";
    file << setw(12) << AcceptanceRatio << "\t";
    file << setw(12) << T << "\n";
    file.close();

}

double Box::computeEnergy() {
    double U = 0;
    for (int i = 0; i < N-1; i++)  {             // all distinct pairs
        for (int j = i+1; j < N; j++) {         // of particles i,j
            double rSqd = MinImage(i, j);
            if (rSqd < rCutOff*rCutOff){
                U += (pow(1 / rSqd, 6) - pow(1 / rSqd, 3));
            }
        }
    }
    return 4 * U;
}

void Box::writeHeaders(ofstream &f, int index, string *fnames) {
    f.open(fnames[index], ios::trunc);
    f.precision(8);
    f << "MCStep s " << setw(12) << "\t " <<
    "Energy_U" << "\t" <<
    setw(12) << "Average_Energy_<U>" << "\t  " <<
    setw(12) << "Variance" << "\t " <<
    setw(12) << "Order_Parameter" << "\t " <<
    setw(12) << "Acceptance_Ratio" << "\t " <<
    setw(12) << "Temperature\n";
}

double *Box::Move(int index, double dr) {
    static double *delta;
    double t = ran0(&l);
    t = 2 * t - 1;
    t = acos(t);
    double ph = ran0(&l);

    delta[0] = dr * cos(t) * sin(2 * ph * PI);
    delta[1] = dr * sin(t) * sin(2 * ph * PI);
    delta[2] = dr * cos(2 * ph * PI);
    Lattice[index][0] += delta[0];
    Lattice[index][1] += delta[1];
    Lattice[index][2] += delta[2];

    return delta;
}

void Box::UnMove(int index, double *delta) {
    Lattice[index][0] -= delta[0];
    Lattice[index][1] -= delta[1];
    Lattice[index][2] -= delta[2];
}

Box::~Box(){
    delete [] Box::Lattice;
    delete [] Box::UAvg;
    delete [] Box::Var;
}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
double ran0(long *idum) {
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
