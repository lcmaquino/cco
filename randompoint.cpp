#include "randompoint.h"

//Note: compile with std=c++11 or std=gnu++11.

RandomPoint::RandomPoint(){
    random_device rd;
    gen.seed(rd());
    min = 0.0;
    max = 1.0;
}

void RandomPoint::generate(double &x, double &y, double &z){
    double t;
    t = dis(gen);
    x = (1.0 - t)*min + t*max;
    t = dis(gen);
    y = (1.0 - t)*min + t*max;
    t = dis(gen);
    z = (1.0 - t)*min + t*max;
}

void RandomPoint::set_limits(double a, double b){
    min = a; max = b;
}
