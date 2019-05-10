#include "dabbene.h"
#include "kharitonov.h"


using namespace std;
using namespace ibex;


void test_topology() {

    Variable x(9);
    IntervalVector x_ini(9);

    double kz     = 0.0;
    double ktheta = 0.0;
    double kpi    = 0.0;
    double kq     = 0.0;
    double iz     = 0.0;
    
    double eps_ctrl = 30;

    double CmQ0  = -0.8116 ;
    double CzW0  = -2.6126 ;
    double CmB10 = -0.7155 ;
    double CzB10 = -1.2233 ;
    double eps   =  0.0000 ;

    int counter = 0;

    x_ini[0] = Interval( -eps_ctrl               , 3*eps_ctrl              );
    x_ini[1] = Interval( -eps_ctrl               , 3*eps_ctrl              );
    x_ini[2] = Interval( -eps_ctrl               , 3*eps_ctrl              );
    x_ini[3] = Interval( -eps_ctrl               , 3*eps_ctrl              );
    x_ini[4] = Interval( -eps_ctrl               , 3*eps_ctrl              );   
    x_ini[5] = Interval( CmQ0  - eps*fabs(CmQ0 ) , CmQ0  + eps*fabs(CmQ0 ) );
    x_ini[6] = Interval( CzW0  - eps*fabs(CzW0 ) , CzW0  + eps*fabs(CzW0 ) );
    x_ini[7] = Interval( CmB10 - eps*fabs(CmB10) , CmB10 + eps*fabs(CmB10) );
    x_ini[8] = Interval( CzB10 - eps*fabs(CzB10) , CzB10 + eps*fabs(CzB10) );

    Function f = Function("Tstab_coefs.txt");
    //Function f = Function(x, -f2(x));

    bool check_stability(Eigen::VectorXd& v);
    
    int iter_test = 100000000;
    IntervalVector  point(9);
    Eigen::VectorXd v(6);

    IntervalVector minmax(9,Interval(0.0,0.0));

    for (int t = 0 ; t < iter_test ; t++) {
        
        for (int j = 0 ; j < 9 ; j++) {
            double val = x_ini.random()[j];
            point[j] = Interval(val, val);
        }
        
        //cout << point << endl;

        IntervalVector res = f.eval_vector(point);
        for (int j = 0 ; j < 6 ; j++) {
            v[j] = res[j].lb();
        }

        //cout << v << endl;
        if(check_stability(v)) {
            counter ++;
            for (int j = 0 ; j < 5 ; j++) {
                if (minmax[j].lb() > point[j].lb()) minmax[j] = Interval(point[j].lb(), minmax[j].ub());
                if (minmax[j].ub() < point[j].lb()) minmax[j] = Interval(minmax[j].lb(), point[j].lb());
            }
        }
    }

    cout << "stables = " << counter << " ; ratio = " << (double) counter / ( (double) iter_test) << endl;
    cout << minmax << endl;
}

int main() {
    test_topology();
    return 1;
}
