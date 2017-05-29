#include "arterialtree.h"
#include "randompoint.h"

class CCO {
private:
    ArterialTree *atree;
    RandomPoint rpoint;
    int N_term, N_con;
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
