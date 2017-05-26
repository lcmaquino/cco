/*
 * arterialtree.h
 *
 *  Created on: 24/05/2017
 *      Author: lcmaquino
 */

#ifndef ARTERIALTREE_H_
#define ARTERIALTREE_H_

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

struct segment{
  int id, up, left, right;
  double beta_l, beta_r, Q, x, y, z, reduced_resistance;
};

class ArterialTree {
private:
    const double viscosity_of_blood = 3.6e-3;
    const double poiseuille_law_constant = 8.0*viscosity_of_blood/M_PI;
    vector<segment> tree;
    int N_term;
    double ox, oy, oz, perfusion_pressure, terminal_pressure,
			perfusion_flow, gamma;
public:
    const int ROOT = 0;
    const int TERMINAL_ENDS = -1;
    ArterialTree(string filename);
    ArterialTree(int Nt, double x, double y, double z, double pperf, double pterm, double Qperf, double gam);
    ~ArterialTree();
    void display();
    void insert(int id, double x, double y, double z);
    void remove(void);
    double flow_splitting_ratio(int id);
    int number_of_terminals(int id);
    double get_radius(int id);
    double get_length(int id);
    bool is_terminal(int id);
    void update(int id);
    vector<int> vicinity(double x, double y, double z, int N_con);
    void save(string filename = "arterialtree.txt");
    void saveVTK(string filename = "arterialtree.vtk");
    void open(string filename);
    void set_origin(double x, double y, double z);
    void insert_root(double x, double y, double z);
    int get_number_of_terminals(void);
    int get_tree_size(void);
    segment get_segment(int id);
    vector<double> get_segement_distal_end(int id);
};

#endif /* ARTERIALTREE_H_ */
