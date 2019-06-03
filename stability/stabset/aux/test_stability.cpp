#include "dabbene.h"
#include "kharitonov.h"
#include "ibex_Random.h"
#include <ctime>

using namespace std;


void test_stability() {
    int iter_test = 1000000;
    int max_iter = 4;
    int degree   = 5;

    ibex::IntervalVector v(1);
    v[0] = ibex::Interval(0,1);

    ibex::IntervalMatrix poly_vect(iter_test,degree+1);
    for (int i = 0 ; i < iter_test ; i ++)
    {
        for (int j = 0 ; j < degree+1 ; j++)
        {   
            double l1 = v.random()[0];
            double l2 = v.random()[0];
            if (l1 < l2) poly_vect[i][j] = ibex::Interval(l1,l1);
            else         poly_vect[i][j] = ibex::Interval(l1,l1);
        }
    }


    for (int t = 0 ; t < iter_test ; t++)
    {
        Eigen::VectorXd poly_p(degree+1);
    
        for (int i = 0 ; i < degree+1 ; i++)
        {
            poly_p(i) = poly_vect[t][i].ub();
        }
    
        bool res_eigen = check_stability_eigen(poly_p);
        bool res_kharitonov = kharitonov(poly_vect.row(t));
        assert(res_eigen == res_kharitonov);
        if (t%10000 == 0) cout << t << endl;
        bool res_dabbene = dabbene(poly_vect.row(t), max_iter);
        if(res_dabbene + res_eigen == 1) cout << " " << res_dabbene << " " << res_eigen << " " << t << endl;
            //if(res_kharitonov) cout << poly << " " << res_kharitonov << " " << res_dabbene << endl;
            //cout << poly << " " << res_kharitonov << " " << res_dabbene << endl;
    }

    //int val = 362465;
    //bool res_kharitonov = kharitonov(poly_vect.row(val));
    //bool res_dabbene = dabbene(poly_vect.row(val), 4);
    //std::cout << res_dabbene << " " << res_kharitonov  << std::endl;
    //std::cout << poly_vect.row(val) << std::endl;
}

int main() {
    ibex::RNG::srand(1234);
    test_stability();
    return 1;
}
