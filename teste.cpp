#include "cco.h"

int main(){
    double pperf = 1.33e4, pterm = 7.98e3, Qperf = 8.33e-6, gamma = 3.0;
    int N_term = 500, N_con = 30;

    CCO c(N_term, N_con, pperf, pterm, Qperf, gamma);
    c.generate_tree();
    c.saveVTK("cco_tree_example.vtk");

    return 0;
}
