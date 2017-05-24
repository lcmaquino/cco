#include "arterialtree.h"
#include <random>

class CCO {
private:
    mt19937 gen;
    uniform_real_distribution<double> dis;
    int N_term, N_con;
    ArterialTree *atree;
    double volume = 1.0, d_min = sqrt(3.0);
public:
    CCO(string filename);
    CCO(int Nt, int Nc, double pperf, double pterm, double Qperf, double gam);
    ~CCO();
    void display();
    void generate_tree(void);
    double evaluate_target_function(void);
    void save(string filename = "cco_tree.txt");
    void saveVTK(string filename = "cco_tree.vtk");
    void open(string filename = "cco_tree.txt");
};
