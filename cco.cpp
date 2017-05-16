#include "cco.h"
#include <fstream>
#include <sstream>
#include <string>

CCO::CCO(int Nt, double pperf, double pterm, double Qperf, double gam) {
    segment root;

    random_device rd;
    gen.seed(rd());

    N_term = Nt;
    tree.reserve(2*N_term - 1);

    ox = dis(gen);
    oy = dis(gen);
    oz = dis(gen);
    gamma = gam;
    perfusion_pressure = pperf;
    terminal_pressure = pterm;

    root.id = 0;
    root.up = -1;
    root.left = -1;
    root.right = -1;
    root.beta_l = 1.0;
    root.beta_r = 1.0;
    root.x = dis(gen);
    root.y = dis(gen);
    root.z = dis(gen);
    root.reduced_resistance = poiseuille_law_constant*sqrt((root.x - ox)*(root.x - ox)
    						  + (root.y - oy)*(root.y - oy) + (root.z - oz)*(root.z - oz));
    root.Q = Qperf;
    tree.push_back(root);

}

CCO::~CCO() {
}

void CCO::insert(int id, double x, double y, double z){
	//note: always insert two elements on the end of tree.
    segment inew, icon;
    double length;

    icon.id = tree.size();
    if (tree[id].left != -1) tree[tree[id].left].up = icon.id;
    if (tree[id].right != -1) tree[tree[id].right].up = icon.id;
    icon.up = id;
    icon.left = tree[id].left;
    icon.right = tree[id].right;
    icon.beta_l = tree[id].beta_l;
    icon.beta_r = tree[id].beta_r;
    icon.Q = tree[id].Q;
    icon.x = tree[id].x;
    icon.y = tree[id].y;
    icon.z = tree[id].z;
    if (tree[id].up != -1) {
    	length = sqrt((tree[id].x - tree[tree[id].up].x)*(tree[id].x - tree[tree[id].up].x)
				  + (tree[id].y - tree[tree[id].up].y)*(tree[id].y - tree[tree[id].up].y) + (tree[id].z - tree[tree[id].up].z)*(tree[id].z - tree[tree[id].up].z));
    }else{
    	length = sqrt((tree[id].x - ox)*(tree[id].x - ox)
				  + (tree[id].y - oy)*(tree[id].y - oy) + (tree[id].z - oz)*(tree[id].z - oz));
    }
    icon.reduced_resistance = tree[id].reduced_resistance - 0.5*poiseuille_law_constant*length;
    tree.push_back(icon);

    tree[id].x /= 2.0;
    tree[id].y /= 2.0;
    tree[id].z /= 2.0;

    inew.id = tree.size();
    inew.up = id;
    inew.left = -1;
    inew.right = -1;
    inew.beta_l = 1.0;
    inew.beta_r = 1.0;
    inew.x = x;
    inew.y = y;
    inew.z = z;
    inew.Q = 1.0;
    inew.reduced_resistance = poiseuille_law_constant*sqrt((x - tree[id].x)*(x - tree[id].x)
    						  + (y - tree[id].y)*(y - tree[id].y) + (z - tree[id].z)*(z - tree[id].z));
    tree.push_back(inew);

    tree[id].left = icon.id;
    tree[id].right = inew.id;

    update(id);
}

void CCO::remove(void){
	int id = tree.size() - 2, id_left = tree[id].left, id_right = tree[id].right,
		id_up = tree[id].up;

	tree[id_up].x *= 2.0;
	tree[id_up].y *= 2.0;
	tree[id_up].z *= 2.0;
	tree[id_up].left = id_left;
	tree[id_up].right = id_right;
	tree[id_left].up = id_up;
	tree[id_right].up = id_up;

	tree.pop_back();
	tree.pop_back();

	update(id_up);
}

int CCO::number_of_terminals(int id){
	if (tree[id].left == -1 && tree[id].right == -1){
		return 1;
	}else{
		return number_of_terminals(tree[id].left) + number_of_terminals(tree[id].right);
	}
}

double CCO::flow_splitting_ratio(int id){
	return (double) number_of_terminals(tree[id].left)/number_of_terminals(tree[id].right);
}

void CCO::update(int id){
	double radii_ratio, radii_ratio_pow, length;
	int i, i_left, i_right, i_up;
    for(i = id; i != -1; i = tree[i].up){
    	i_left = tree[i].left;
    	i_right = tree[i].right;
    	i_up = tree[i].up;
        if (i_up != -1) {
        	length = sqrt((tree[i].x - tree[i_up].x)*(tree[i].x - tree[i_up].x)
    				  + (tree[i].y - tree[i_up].y)*(tree[i].y - tree[i_up].y)
					  + (tree[i].z - tree[i_up].z)*(tree[i].z - tree[i_up].z));
        }else{
        	length = sqrt((tree[i].x - ox)*(tree[i].x - ox)
    				  + (tree[i].y - oy)*(tree[i].y - oy)
					  + (tree[i].z - oz)*(tree[i].z - oz));
        }
    	if (i_left == -1 && i_right == -1){
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

void CCO::generate_tree(void){
	double x, y, z, smin_dist, sdist, tol = 0.0000000001;
	int id_min_dist;
	for (int i = 0; i < N_term - 1; i++){
		x = dis(gen);
		y = dis(gen);
		z = dis(gen);
		smin_dist = 1.0;
		id_min_dist = 0;
		for(vector<segment>::iterator it = tree.begin(); it != tree.end(); ++it){
			sdist = ((*it).x - x)*((*it).x - x) + ((*it).y - y)*((*it).y - y)
					+ ((*it).z - z)*((*it).z - z);
			if (fabs(sdist - smin_dist) < tol) {
				id_min_dist = (*it).id;
				smin_dist = sdist;
			}
		}
		insert(id_min_dist, x, y, z);
	}
}

void CCO::display(){
    for(vector<segment>::iterator it = tree.begin(); it != tree.end(); ++it){
        cout << "(id = " << (*it).id << ", (" << (*it).x <<  ", " << (*it).y <<  ", "  
             << (*it).z << "), " << (*it).up << ", " << (*it).left << ", "
			 << (*it).right << ", " << (*it).reduced_resistance << ", "
			 << (*it).beta_l << ", " << (*it).beta_r << " )" << endl;
    }
}

void CCO::save(string filename){
    ofstream treefile;
    treefile.open(filename);
    if (treefile.is_open()) {
        treefile << N_term << " " << ox << " " << oy << " " << oz << " "
        << perfusion_pressure << " " << terminal_pressure << " "
        << gamma << endl;
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

void CCO::open(string filename){
    ifstream treefile;
    string line;
    treefile.open(filename);
    if (treefile.is_open()) {
        getline(treefile, line);
        istringstream iss(line);
        iss >> N_term >> ox >> oy >>  oz >> perfusion_pressure
        >> terminal_pressure >> gamma;
        int i = 0;
        while (getline(treefile, line)) {
            istringstream iss(line);
            if (i == 0) {
                iss >> tree[0].x >> tree[0].y >> tree[0].z >> tree[0].up >> tree[0].left >> tree[0].right >> tree[0].reduced_resistance >> tree[0].beta_l >> tree[0].beta_r;
            }else{
                segment s;
                s.id = i;
                iss >> s.x >> s.y >> s.z >> s.up >> s.left >> s.right >> s.reduced_resistance >> s.beta_l >> s.beta_r;
                tree.push_back(s);
            }
            i++;
        }
        treefile.close();
    }else{
        cout << "Unable to open \"" << filename << "\"." << endl;
    }
}
