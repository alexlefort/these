#include "optim_KD.h"
#include "ibex_Random.h"
#include "ibex.h"
#include <ctime>

using namespace std;
using namespace ibex;

void test_optim_kd() {

	Variable x(5);
    IntervalVector x_ini(5);
    
    double eps_ctrl  = 10.0     ; 
    double eps       =  0.00001 ;
    double CzW0      = -2.7     ;
    double CmQ0      = -0.38    ;
    
    x_ini[0] = Interval(-eps_ctrl,eps_ctrl);
    x_ini[1] = Interval(-eps_ctrl,eps_ctrl);
    x_ini[2] = Interval(-eps_ctrl,eps_ctrl);
    x_ini[3] = Interval(-eps_ctrl,eps_ctrl);
    x_ini[4] = Interval(-eps_ctrl,eps_ctrl);
    //x_ini[5] = Interval(CzW0 - eps*fabs(CzW0), CzW0 + eps*fabs(CzW0));
    //x_ini[6] = Interval(CmQ0 - eps*fabs(CmQ0), CmQ0 + eps*fabs(CmQ0));


    Function coeffs = Function("Tstab_coefs.txt");

    const ExprNode& stab1 = coeffs(x)[0];
    const ExprNode& stab2 = coeffs(x)[1];
    const ExprNode& stab3 = coeffs(x)[2];
    const ExprNode& stab4 = coeffs(x)[3];
    const ExprNode& stab5 = coeffs(x)[4];
    const ExprNode& stab6 = coeffs(x)[5];

    Function stab(x,ibex::max(ibex::max(ibex::max(ibex::max(-stab6,-stab4),-stab2),ibex::pow(stab1,2)*ibex::pow(stab6,2) + ibex::pow(stab2,2)*ibex::pow(stab5,2) + stab1*ibex::pow(stab4,2)*stab5 + stab2*ibex::pow(stab3,2)*stab6 - 2*stab1*stab2*stab5*stab6 - stab1*stab3*stab4*stab6 - stab2*stab3*stab4*stab5),stab1*stab4 - stab2*stab3)); 

    SystemFactory fac_x;
    fac_x.add_var(x, x_ini);
    fac_x.add_goal(coeffs);

    NormalizedSystem sys(fac_x); 

    int nb_iter = 1000000;
    int trace_freq = 1000;
    cout << "start" << endl;
    OptimKD oo(sys, x_ini, nb_iter, trace_freq);

    oo.optimize();

}

int main() {
    ibex::RNG::srand(1234);
    test_optim_kd();
    return 1;
}