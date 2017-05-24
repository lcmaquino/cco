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

	//do {
		x = dis(gen);
		y = dis(gen);
		z = dis(gen);
	//}while(isnan(x) || isnan(y) || isnan(z));

    //cout << "Origin: " << x << ", " << y << ", " << z << endl;
	atree = new ArterialTree(Nt, x, y, z, pperf, pterm, Qperf, gam);

	//do {
		x = dis(gen);
		y = dis(gen);
		z = dis(gen);
	//}while(isnan(x) || isnan(y) || isnan(z));

    //cout << "First root: " << x << ", " << y << ", " << z << endl;

	atree->insert_root(x, y, z);
}

CCO::~CCO() {
	delete atree;
}

void CCO::generate_tree(void){
	double x, y, z;
	vector<int> vic;
	int id, k, ntry = 20;
	double T, T_min, dtemp;
	for (int i = 0; i < N_term - 1; i++){
		k = 0;
		do {
			vector<double> dist;
			x = dis(gen);
			y = dis(gen);
			z = dis(gen);
		//cout << "Random Point: " << x << ", " << y << ", " << z << endl;
			vic = atree->vicinity(x, y, z, N_con, dist);
			dtemp = dist[0];
			/*
			if (k == 0) {
				for (int v : vic) {
					cout << v << " ";
				}
				cout << endl;
				for (double d : dist) {
					cout << d << " ";
				}
				cout << endl;
			}*/
			k++;
		}while(dtemp < d_min && k < ntry);
		//cout << "Minimum distance: " << dtemp << " d_min = " << d_min << " Tries: " << k << endl;
		//}while(k < ntry);
		for(unsigned int j = 0; j < vic.size(); j++){
			//cout << vic[j] << " ";
			if (j == 0){
				id = vic[0];
				atree->insert(vic[0], x, y, z);
				T_min = evaluate_target_function();
				//cout << T_min << " ";
			}else{
				atree->insert(vic[j], x, y, z);
				T = evaluate_target_function();
				if (T < T_min) {
					id = vic[j];
					T = T_min;
				}
				//cout << T << " ";
			}
			atree->remove();
		}
		//cout << endl;
		//cout << "id = " << id << " T_min = " << T_min << endl;
		atree->insert(id, x, y, z);
		d_min = sqrt(3.0)*pow(volume/atree->number_of_terminals(atree->ROOT), 1.0/3.0);
		//atree->display();
		//cout << "---" << endl;
	}
}

double CCO::evaluate_target_function(void){
	int i, size = atree->get_tree_size();
	double T = 0.0, radius, length;
	for(i = atree->ROOT; i < size; i++){
		radius = atree->get_radius(i);
		length = atree->get_length(i);
		//if (isnan(radius) || isnan(length))
		//	cout << i << ", " <<radius << ", " << length << endl;
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
