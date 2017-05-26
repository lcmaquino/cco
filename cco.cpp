#include "cco.h"
#include <fstream>
#include <sstream>
#include <string>

CCO::CCO(string filename) {
	atree = new ArterialTree(filename);
	N_term = atree->get_number_of_terminals();
	N_con = 10;
}

CCO::CCO(int Nt, int Nc, double pperf, double pterm, double Qperf, double gam) {
	double x, y, z;
	random_device rd;

	N_term = Nt;
	N_con = Nc;
	gen.seed(rd());
	//Generate a random origin.
	x = 0.0; //x = dis(gen);
	y = 0.0; //y = dis(gen);
	z = 0.0; //z = dis(gen);
	atree = new ArterialTree(Nt, x, y, z, pperf, pterm, Qperf, gam);
	//Generate a random root distal end.
	x = dis(gen);
	y = dis(gen);
	z = dis(gen);
	atree->insert_root(x, y, z);
}

CCO::~CCO() {
	delete atree;
}

void CCO::generate_tree(void){
	double x, y, z;
	vector<int> vic;
	vector<double> point;
	int id, att, attempts = 20;
	double d, T_min, T;
	for (int i = 0; i < N_term - 1; i++){
		att = 0;
		//Try to generate a random new segment that has
		//at least d_min distance of all segments generated
		//so far.
		do {
			x = dis(gen);
			y = dis(gen);
			z = dis(gen);
			//Get the N_con segments near to (x, y ,z).
			vic = atree->vicinity(x, y, z, N_con);
			//vic[0] is nearest segment.
			point = atree->get_segement_distal_end(vic[0]);
			d = sqrt( (x - point[0])*(x - point[0])
					+ (y - point[1])*(y - point[1])
					+ (z - point[2])*(z - point[2]));
			att++;
		}while(d < d_min && att < attempts);

		//Find where to connect the new segment
		//so the target function is minimum.
		id = vic[0];
		atree->insert(vic[0], x, y, z);
		T_min = evaluate_target_function();
		atree->remove();
		for(unsigned int j = 1; j < vic.size(); j++){
			atree->insert(vic[j], x, y, z);
			T = evaluate_target_function();
			if (T < T_min) {
				id = vic[j];
				T = T_min;
			}
			atree->remove();
		}

		//Insert the new segment.
		atree->insert(id, x, y, z);
		//Update d_min.
		d_min = sqrt(3.0)*pow(volume/atree->number_of_terminals(atree->ROOT), 1.0/3.0);
	}
}

double CCO::evaluate_target_function(void){
	unsigned int i, size = atree->get_tree_size();
	double T = 0.0, radius, length;
	for(i = atree->ROOT; i < size; i++){
		radius = atree->get_radius(i);
		length = atree->get_length(i);
		T += length*radius*radius;
	}
	return T;
}

void CCO::display(){
	atree->display();
}

void CCO::saveVTK(string filename){
	atree->saveVTK(filename);
}

void CCO::save(string filename){
	atree->save(filename);
}

void CCO::open(string filename){
	atree->open(filename);
}
