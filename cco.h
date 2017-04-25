#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

struct segment{
  int id, left, right;
  double beta_l, beta_r, Q, x, y, z;
};

class CCO {
private:
    vector<segment> tree;
    int segment_count;
public:
    CCO(double Q, double x, double y, double z);
    ~CCO();
    void display();
    void insert(int id, double x, double y, double z);
};
