/*
 * arterialtree.cpp
 *
 *  Created on: 24/05/2017
 *      Author: lcmaquino
 */

#include "arterialtree.h"
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

ArterialTree::ArterialTree(string filename) {
	open(filename);
}

ArterialTree::ArterialTree(int Nt, double x, double y, double z,
	     double pperf, double pterm, double Qperf, double gam){
    ox = x; oy = y; oz = z;
    N_term = Nt;
    tree.reserve(2*N_term - 1);
    gamma = gam;
    perfusion_pressure = pperf;
    terminal_pressure = pterm;
    perfusion_flow = Qperf;
}

ArterialTree::~ArterialTree() {
}

void ArterialTree::insert(int id, double x, double y, double z){
	//note: always insert two elements on the end of tree.
    segment inew, icon;
    double length;

    icon.id = tree.size();
    if (tree[id].left != TERMINAL_ENDS) tree[tree[id].left].up = icon.id;
    if (tree[id].right != TERMINAL_ENDS) tree[tree[id].right].up = icon.id;
    icon.up = id;
    icon.left = tree[id].left;
    icon.right = tree[id].right;
    icon.beta_l = tree[id].beta_l;
    icon.beta_r = tree[id].beta_r;
    icon.x = tree[id].x;
    icon.y = tree[id].y;
    icon.z = tree[id].z;
    length = get_length(id);
    icon.reduced_resistance = tree[id].reduced_resistance - 0.5*poiseuille_law_constant*length;
    tree.push_back(icon);

    if (id == ROOT) {
        tree[id].x = (tree[id].x + ox)/2.0;
        tree[id].y = (tree[id].y + oy)/2.0;
        tree[id].z = (tree[id].z + oz)/2.0;
    }else{
        tree[id].x = (tree[id].x + tree[tree[id].up].x)/2.0;
        tree[id].y = (tree[id].y + tree[tree[id].up].y)/2.0;
        tree[id].z = (tree[id].z + tree[tree[id].up].z)/2.0;
    }

    inew.id = tree.size();
    inew.up = id;
    inew.left = TERMINAL_ENDS;
    inew.right = TERMINAL_ENDS;
    inew.beta_l = 1.0;
    inew.beta_r = 1.0;
    inew.x = x;
    inew.y = y;
    inew.z = z;
    inew.reduced_resistance = poiseuille_law_constant*sqrt((x - tree[id].x)*(x - tree[id].x)
    						  + (y - tree[id].y)*(y - tree[id].y) + (z - tree[id].z)*(z - tree[id].z));
    tree.push_back(inew);

    tree[id].left = icon.id;
    tree[id].right = inew.id;

    update(id);
}

void ArterialTree::remove(void){
	int id = tree.size() - 2, id_left = tree[id].left, id_right = tree[id].right,
		id_up = tree[id].up;
	/*
    if (id == 1) {
    	tree[id_up].x *= 2.0;
    	tree[id_up].y *= 2.0;
    	tree[id_up].z *= 2.0;
    	tree[id_up].x -= ox;
    	tree[id_up].y -= oy;
    	tree[id_up].z -= oz;
    }else{
    	tree[id_up].x *= 2.0;
    	tree[id_up].y *= 2.0;
    	tree[id_up].z *= 2.0;
    	tree[id_up].x -= tree[tree[id_up].up].x;
    	tree[id_up].y -= tree[tree[id_up].up].y;
    	tree[id_up].z -= tree[tree[id_up].up].z;
    }
    */

	tree[id_up].x = tree[id].x;
	tree[id_up].y = tree[id].y;
	tree[id_up].z = tree[id].z;

	tree[id_up].left = id_left;
	tree[id_up].right = id_right;

	if (id_left != -1) tree[id_left].up = id_up;
	if (id_right != -1) tree[id_right].up = id_up;

	tree.pop_back();
	tree.pop_back();

	update(id_up);
}

bool ArterialTree::is_terminal(int id){
	return (tree[id].left == TERMINAL_ENDS && tree[id].right == TERMINAL_ENDS);
}

int ArterialTree::number_of_terminals(int id){
	if (is_terminal(id)){
		return 1;
	}else{
		return number_of_terminals(tree[id].left) + number_of_terminals(tree[id].right);
	}
}

double ArterialTree::flow_splitting_ratio(int id){
	return (double) number_of_terminals(tree[id].left)/number_of_terminals(tree[id].right);
}

void ArterialTree::update(int id){
	double radii_ratio, radii_ratio_pow, length;
	int i, i_left, i_right;
    for(i = id; i != TERMINAL_ENDS; i = tree[i].up){
    	i_left = tree[i].left;
    	i_right = tree[i].right;
        length = get_length(i);
    	if (is_terminal(i)){
    		tree[i].beta_l = 1.0;
    		tree[i].beta_r = 1.0;
            tree[i].reduced_resistance = poiseuille_law_constant*length;
    	}else{
        	radii_ratio = pow(flow_splitting_ratio(i)*tree[i_left].reduced_resistance/tree[i_right].reduced_resistance, 0.25);
            radii_ratio_pow = pow(radii_ratio, gamma);
            tree[i].beta_l = pow(1 + 1/radii_ratio_pow, -1/gamma);
            tree[i].beta_r = pow(1 + radii_ratio_pow, -1/gamma);
            tree[i].reduced_resistance = poiseuille_law_constant*length
            							  + 1/(pow(tree[i].beta_l, 4.0)/tree[i_left].reduced_resistance
        								  + pow(tree[i].beta_r, 4.0)/tree[i_right].reduced_resistance);
    	}
    }
}

vector<int> ArterialTree::vicinity(double x, double y, double z, int N_con){
	struct distance_index{
		int id;
		double dist;
	};

    struct {
        bool operator()(distance_index a, distance_index b) const
        {
            return a.dist < b.dist;
        }
    }comparison_function;

	distance_index aux;
	vector<distance_index> di;
	vector<int> v;
	int n = N_con < tree.size() ? N_con : tree.size();

	for(vector<segment>::iterator it = tree.begin(); it != tree.end(); ++it){
		aux.id = (*it).id;
		//We just need sort the distance in ascendent order.
		//We can sort the values of (x - x0)^2 + (y - y0)^2 + (z - z0)^2
		//instead of sqrt((x - x0)^2 + (y - y0)^2 + (z - z0)^2).
		//This way we use one less operation.
		aux.dist = ((*it).x - x)*((*it).x - x) + ((*it).y - y)*((*it).y - y)
				+ ((*it).z - z)*((*it).z - z);
		di.push_back(aux);
	}

	partial_sort(di.begin(), di.begin() + n, di.end(), comparison_function);

	for(int i = 0; i < n; i++){
		v.push_back(di[i].id);
	}
	return v;
}

void ArterialTree::display(){
    for(vector<segment>::iterator it = tree.begin(); it != tree.end(); ++it){
        cout << "(id = " << (*it).id << ", (" << (*it).x <<  ", " << (*it).y <<  ", "
             << (*it).z << "), " << (*it).up << ", " << (*it).left << ", "
			 << (*it).right << ", " << (*it).reduced_resistance << ", "
			 << (*it).beta_l << ", " << (*it).beta_r << " )" << endl;
    }
}

void ArterialTree::set_origin(double x, double y, double z){
	ox = x;
	oy = y;
	oz = z;
}

void ArterialTree::insert_root(double x, double y, double z){
	segment root;
    root.id = ROOT;
    root.up = TERMINAL_ENDS;
    root.left = TERMINAL_ENDS;
    root.right = TERMINAL_ENDS;
    root.x = x; root.y = y; root.z = z;
    tree.push_back(root);
    update(ROOT);
}

double ArterialTree::get_length(int id){
	double length;
	int id_up;
	id_up = tree[id].up;
	if (id_up != TERMINAL_ENDS) {
		length = sqrt((tree[id].x - tree[id_up].x)*(tree[id].x - tree[id_up].x)
				  + (tree[id].y - tree[id_up].y)*(tree[id].y - tree[id_up].y)
				  + (tree[id].z - tree[id_up].z)*(tree[id].z - tree[id_up].z));
	}else{
		length = sqrt((tree[id].x - ox)*(tree[id].x - ox)
				  + (tree[id].y - oy)*(tree[id].y - oy)
				  + (tree[id].z - oz)*(tree[id].z - oz));
	}
	return length;
}

double ArterialTree::get_radius(int id){
    double r_root;
    if (id == ROOT) {
    	r_root = pow(tree[ROOT].reduced_resistance*perfusion_flow/
    	    					(perfusion_pressure - terminal_pressure), 0.25);
    	return r_root;
    }else{
    	if (tree[tree[id].up].left == id)
    		return tree[tree[id].up].beta_l*get_radius(tree[id].up);
    	else
    		return tree[tree[id].up].beta_r*get_radius(tree[id].up);
    }
}

void ArterialTree::saveVTK(string filename){
    ofstream treefile;
    int N_segments = tree.size();
    treefile.open(filename);
    if (treefile.is_open()) {
        treefile << "# vtk DataFile Version 3.0" << endl;
        treefile << "vtk output" << endl;
        treefile << "ASCII" << endl;
        treefile << "DATASET POLYDATA" << endl;
        treefile << "POINTS  " << N_segments + 1 << "  float" << endl;
        treefile << ox << "  " << oy << "  " << oz << endl;
        for(vector<segment>::iterator it = tree.begin(); it != tree.end(); ++it){
            treefile << (*it).x << "  " << (*it).y << "  " << (*it).z << endl;
        }
        treefile << endl;
        treefile << "LINES  " << N_segments << "  " << 3*N_segments << endl;
        treefile << "2  0  1" << endl;
        for(vector<segment>::iterator it = tree.begin(); it != tree.end(); ++it){
        	if (!is_terminal((*it).id)) {
                treefile << "2  " << (*it).id + 1 << "  " << (*it).left + 1 << endl;
                treefile << "2  " << (*it).id + 1 << "  " << (*it).right + 1 << endl;
        	}
        }
        treefile << endl;
		treefile << "CELL_DATA  " << N_segments << endl;
		treefile << "scalars radius float" << endl;
		treefile << "LOOKUP_TABLE default" << endl;
		treefile << get_radius(ROOT) << endl;
        for(vector<segment>::iterator it = tree.begin(); it != tree.end(); ++it){
        	if (!is_terminal((*it).id)) {
                treefile << get_radius((*it).left) << endl;
                treefile << get_radius((*it).right) << endl;
        	}
        }
        treefile.close();
    }else{
        cout << "Unable to open \"" << filename << "\"." << endl;
    }
}

void ArterialTree::save(string filename){
    ofstream treefile;
    treefile.open(filename);
    if (treefile.is_open()) {
        treefile << N_term << " "
        		 << ox << " " << oy << " " << oz << " "
				 << perfusion_pressure << " " << terminal_pressure << " "
				 << perfusion_flow << " " << gamma << endl;
        for(vector<segment>::iterator it = tree.begin(); it != tree.end(); ++it){
            treefile << (*it).x << " " << (*it).y << " " << (*it).z << " "
            << (*it).up << " " << (*it).left << " " << (*it).right << " "
            << (*it).reduced_resistance  << " " << (*it).beta_l << " "
            << (*it).beta_r  << " " << endl;
        }
        treefile.close();
    }else{
        cout << "Unable to open \"" << filename << "\"." << endl;
    }
}

void ArterialTree::open(string filename){
    ifstream treefile;
    string line;
    treefile.open(filename);
    if (treefile.is_open()) {
        getline(treefile, line);
        istringstream iss(line);
        iss >> N_term >> ox >> oy >>  oz >> perfusion_pressure
        >> terminal_pressure >> perfusion_flow >> gamma;
        tree.reserve(2*N_term - 1);
        int i = 0;
        while (getline(treefile, line)) {
            istringstream iss(line);
			segment s;
			s.id = i;
			iss >> s.x >> s.y >> s.z >> s.up >> s.left >> s.right >> s.reduced_resistance >> s.beta_l >> s.beta_r;
			tree.push_back(s);
            i++;
        }
        treefile.close();
    }else{
        cout << "Unable to open \"" << filename << "\"." << endl;
    }
}

int ArterialTree::get_number_of_terminals(void){
	return N_term;
}

int ArterialTree::get_tree_size(void){
	return tree.size();
}

segment ArterialTree::get_segment(int id){
	return tree[id];
}

vector<double> ArterialTree::get_segement_distal_end(int id){
	vector<double> v = {tree[id].x, tree[id].y, tree[id].z};
	return v;
}


