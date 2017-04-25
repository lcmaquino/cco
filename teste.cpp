#include "cco.h"

int main(){
    double pperf = 100.0, pterm = 60.0, Qperf = 500.0, gamma = 3.0;

    CCO c(Qperf, 1.0, 1.0, 1.0);

    c.display();
    cout << "---" << endl;

    c.insert(0, 2.2, 2.4, 2.6);
    c.display();
    cout << "---" << endl;;

    c.insert(0, 3.2, 3.4, 3.6);
    c.display();
    cout << "---" << endl;

    c.insert(2, 4.2, 4.4, 4.6);
    c.display();
    cout << "---" << endl;

    return 0;
}
