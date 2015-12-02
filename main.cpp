#include "Box.h"

using namespace std;
const int MCSteps = 100;
const int N = 32;               // number of particles
int main() {
    time_t s_i;
    s_i = time(NULL);
    string filenames_A[3] = {"fcc_a0.txt", "fcc_a1.txt", "fcc_a2.txt"};
    string filenames_B[3] = {"fcc_b0.txt", "fcc_b1.txt", "fcc_b2.txt"};
    double temperatures[3] = {0.1, 0.5, 1};
    double stepSizes[3] = {0.06, 0.065, 0.08};
    ofstream file_a;
    ofstream file_b;
    Box *BoxA;
    Box *BoxB;
    for (int i = 0; i < 3; i++) {
        BoxA = new Box(N, MCSteps, "Right");
        BoxB = new Box(N, MCSteps, "Left");
        BoxA->writeHeaders(file_a, i, filenames_A);
        BoxB->writeHeaders(file_b, i, filenames_B);
        for (int s = 0; s < MCSteps; s++) {
            BoxA->MMC(s, temperatures[i], filenames_A[i], stepSizes[i]);
            BoxB->MMC(s, temperatures[i], filenames_B[i], stepSizes[i]);
        }
        delete BoxA;
        delete BoxB;
    }

    time_t s_f;
    s_f = time(NULL);
    cout << "Runtime : " << s_f - s_i <<" seconds \n";
    //system ("pause");
    return 0;
}