#include <random>
#include <iostream>

using namespace std;

class RandomPoint{
    private:
    mt19937 gen;
    uniform_real_distribution<double> dis;
    double min, max;
    public:
    RandomPoint();
    void generate(double &x, double &y, double &z);
    void set_limits(double a, double b);
};

