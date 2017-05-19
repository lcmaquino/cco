#include "cco.h"

int main(){
    double pperf = 100.0, pterm = 60.0, Qperf = 500.0, gamma = 3.0;
    int N_term = 4000;

    CCO c(N_term, pperf, pterm, Qperf, gamma);
    c.generate_tree();
    //CCO c("cco_tree_5_exemplo1.txt");
    c.display();
    cout << "---" << endl;

    c.save("cco_tree_4000_exemplo1.txt");
    c.saveVTK("cco_tree_4000_exemplo1.vtk");

    /*
    c.set_origin(0.0, 1.0, 0.0);
    c.set_root_end(0.0, 0.0, 0.0);
    c.display();
    cout << "---" << endl;

    c.insert(0, 0.5, 0.0, 0.0);
    c.display();
    cout << "---" << endl;

    c.insert(1, -0.4, -0.6, 0.0);
    c.display();
    cout << "---" << endl;

    c.insert(2, 0.4, -0.2, 0.0);
    c.display();
    cout << "---" << endl;

    c.remove();
    c.display();
    cout << "---" << endl;

    c.remove();
    c.display();
    cout << "---" << endl;

    c.remove();
    c.display();
    cout << "---" << endl;
    */
    //c.save("cco_tree_3_exemplo2.txt");

    //CCO c("cco_tree_9_exemplo2.txt");
    //c.saveVTK("cco_tree_3_exemplo2.vtk");


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
