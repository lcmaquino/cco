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
    vector<point> points;
    mt19937 gen;
    uniform_real_distribution<> dis;
    int segment_count, N_term;
    double ox, oy, oz, perfusion_pressure, terminal_pressure, gamma;
    const double viscosity_of_blood = 3.6;
    const double poiseuille_law_constant = 8.0*viscosity_of_blood/M_PI;
    int number_of_terminals(int id);
public:
    CCO(int Nt, double pperf, double pterm, double Qperf, double gam);
    ~CCO();
    void display();
    void insert(int id, double x, double y, double z);
    void remove(void);
    double flow_splitting_ratio(int id);
    void update(int id);
    void generate_tree(void);
};
