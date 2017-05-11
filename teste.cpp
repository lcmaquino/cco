#include "cco.h"

int main(){
    double pperf = 100.0, pterm = 60.0, Qperf = 500.0, gamma = 3.0;
    int N_term = 5;

    CCO c(N_term, pperf, pterm, Qperf, gamma);

    c.generate_tree();
    c.display();
    cout << "---" << endl;
    /*
    c.insert(0, 2.2, 2.4, 2.6);
    c.display();
    cout << "---" << endl;;

    c.remove();
    c.display();
    cout << "---" << endl;;

    c.insert(0, 2.2, 2.4, 2.6);
    c.display();
    cout << "---" << endl;;

    c.insert(0, 3.2, 3.4, 3.6);
    c.display();
    cout << "---" << endl;

    c.insert(2, 4.2, 4.4, 4.6);
    c.display();
    cout << "---" << endl;

    c.insert(3, 5.2, 5.4, 5.6);
    c.display();
    cout << "---" << endl;

    c.remove();
    c.display();
    cout << "---" << endl;

    c.remove();
    c.display();
    cout << "---" << endl;
	*/
    return 0;
}
