#include "Box.h"

using namespace std;
const int MCSteps = 100;
const int N = 32;               // number of particles
int main() {
    time_t s_i;
    s_i = time(NULL);
    string filenames_A[3] = {"fcc_a0.txt", "fcc_a1.txt", "fcc_a2.txt"};
    string filenames_B[3] = {"fcc_b0.txt", "fcc_b1.txt", "fcc_b2.txt"};
    if (remove("fcc_a0.txt") == 0
        && remove("fcc_a1.txt") == 0
        && remove("fcc_a2.txt") == 0
        && remove("fcc_a2.txt") == 0
        && remove("fcc_a2.txt") == 0
        && remove("fcc_a2.txt") == 0
            )
        puts("Files successfully deleted");
    else perror("Error deleting file");
    double temperatures[3] = {0.1, 0.5, 1};
    double steps[3] = {0.06, 0.065, 0.08};
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
            BoxA->MMSteps(s, temperatures[i], filenames_A[i], steps[i]);
            BoxB->MMSteps(s, temperatures[i], filenames_A[i], steps[i]);
        }
        delete BoxA;
        delete BoxB;
    }

    time_t s_f;
    s_f = time(NULL);
    cout << "Runtime : " << s_f - s_i <<" seconds \n";
    system ("pause");
    return 0;
}