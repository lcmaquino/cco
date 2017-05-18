#include <iostream>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

struct point{
	double x, y, z;
};

struct segment{
  int id, up, left, right;
  double beta_l, beta_r, Q, x, y, z, reduced_resistance;
};

class CCO {
private:
    vector<segment> tree;
    mt19937 gen;
    uniform_real_distribution<> dis;
    int N_term;
    double ox, oy, oz, perfusion_pressure, terminal_pressure, perfusion_flow, gamma;
    const double viscosity_of_blood = 3.6;
    const double poiseuille_law_constant = 8.0*viscosity_of_blood/M_PI;
    int number_of_terminals(int id);
    double get_radius(int id);
    double get_length(int id);
public:
    CCO(string filename);
    CCO(int Nt, double pperf, double pterm, double Qperf, double gam);
    ~CCO();
    void display();
    void insert(int id, double x, double y, double z);
    void remove(void);
    double flow_splitting_ratio(int id);
    void update(int id);
    void generate_tree(void);
    void save(string filename = "cco_tree.txt");
    void saveVTK(string filename = "cco_tree.vtk");
    void open(string filename = "cco_tree.txt");
};
