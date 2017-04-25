#include "cco.h"

CCO::CCO(double Q, double x, double y, double z) {
    segment root;
    root.id = 0;
    root.up = -1;
    root.left = -1;
    root.right = -1;
    root.beta_l = 1.0;
    root.beta_r = 1.0;
    root.x = x;
    root.y = y;
    root.z = z;
    root.Q = Q;
    tree.push_back(root);
    segment_count = 1;
}

CCO::~CCO() {
}

void CCO::insert(int id, double x, double y, double z){

    segment inew, ibif, icon;

    icon.id = tree.size();
    icon.up = id;
    icon.left = tree[id].left;
    if (tree[id].left != -1) tree[tree[id].left].id = icon.id;
    icon.right = tree[id].right;
    if (tree[id].right != -1) tree[tree[id].right].id = icon.id;
    icon.beta_l = tree[id].beta_l;
    icon.beta_r = tree[id].beta_r;
    icon.Q = tree[id].Q;
    icon.x = tree[id].x;
    icon.y = tree[id].y;
    icon.z = tree[id].z;
    tree.push_back(icon);

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
    tree.push_back(inew);

    tree[id].left = icon.id;
    tree[id].right = inew.id;
    tree[id].x /= 2.0;
    tree[id].y /= 2.0;
    tree[id].z /= 2.0;

    segment_count += 2;
}

void CCO::display(){
    for(vector<segment>::iterator it = tree.begin(); it != tree.end(); ++it){
        cout << "(id = " << (*it).id << ", (" << (*it).x <<  ", " << (*it).y <<  ", "  
             << (*it).z << "), " << (*it).up << ", " << (*it).left << ", " << (*it).right << " )" << endl;
    }
}
