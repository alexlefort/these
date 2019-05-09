#include "ibex.h"
#include <string>
#include <ctime>
#include "ibex_Random.h"


using namespace ibex;
//using namespace std;

void test_stability() {

    IntervalVector x(9);

    double kz     = 0.0;
    double ktheta = 0.0;
    double kpi    = 0.0;
    double kq     = 0.0;
    double iz     = 0.0;
    
    double eps_ctrl = 5.0;

    x[0] = Interval(-eps_ctrl,   eps_ctrl);
    x[1] = Interval(-eps_ctrl, 3*eps_ctrl);
    x[2] = Interval(-eps_ctrl, 3*eps_ctrl);
    x[3] = Interval(-1       , 3*eps_ctrl);
    x[4] = Interval(-eps_ctrl,   eps_ctrl);

    double CmQ0  = -0.8116 ;
    double CzW0  = -2.6126 ;
    double CmB10 = -0.7155 ;
    double CzB10 = -1.2233 ;
    
    x[5] = Interval(CmQ0 , CmQ0 );
    x[6] = Interval(CzW0 , CzW0 );
    x[7] = Interval(CmB10, CmB10);
    x[8] = Interval(CzB10, CzB10);

    Variable p(9); 

    Function *coeffs = new Function("functions/Tstab_coefs_x.txt");

    const ExprNode& stab1 = (*coeffs)(p)[0];
    const ExprNode& stab2 = (*coeffs)(p)[1];
    const ExprNode& stab3 = (*coeffs)(p)[2];
    const ExprNode& stab4 = (*coeffs)(p)[3];
    const ExprNode& stab5 = (*coeffs)(p)[4];
    const ExprNode& stab6 = (*coeffs)(p)[5];
    const ExprNode& stab7 = (*coeffs)(p)[6];

    Vector x_min(5), x_max(5);

    x_min[0] =  eps_ctrl;
    x_min[1] =  eps_ctrl;
    x_min[2] =  eps_ctrl;
    x_min[3] =  eps_ctrl;
    x_min[4] =  eps_ctrl;

    x_max[0] = -eps_ctrl;
    x_max[1] = -eps_ctrl;
    x_max[2] = -eps_ctrl;
    x_max[3] = -eps_ctrl;
    x_max[4] = -eps_ctrl;

    Function* stab = new Function(p,ibex::max(ibex::max(ibex::max(ibex::max(ibex::max(ibex::max(-stab7,-stab5),-stab3),-stab1),ibex::pow(stab2,3)*ibex::pow(stab7,3) - stab1*ibex::pow(stab4,3)*ibex::pow(stab7,2) + ibex::pow(stab1,2)*ibex::pow(stab6,3)*stab7 + stab2*stab3*ibex::pow(stab4,2)*ibex::pow(stab7,2) + stab2*ibex::pow(stab3,2)*ibex::pow(stab6,2)*stab7 - 2*ibex::pow(stab2,2)*stab3*stab6*ibex::pow(stab7,2) - ibex::pow(stab2,2)*stab4*stab5*ibex::pow(stab7,2) + ibex::pow(stab2,2)*ibex::pow(stab5,2)*stab6*stab7 + 3*stab1*stab2*stab4*stab6*ibex::pow(stab7,2) - 2*stab1*stab2*stab5*ibex::pow(stab6,2)*stab7 - stab1*stab3*stab4*ibex::pow(stab6,2)*stab7 + stab1*ibex::pow(stab4,2)*stab5*stab6*stab7 - stab2*stab3*stab4*stab5*stab6*stab7),ibex::pow(stab1,2)*ibex::pow(stab6,2) + ibex::pow(stab2,2)*ibex::pow(stab5,2) + stab1*ibex::pow(stab4,2)*stab5 + stab2*ibex::pow(stab3,2)*stab6 - ibex::pow(stab2,2)*stab3*stab7 + stab1*stab2*stab4*stab7 - 2*stab1*stab2*stab5*stab6 - stab1*stab3*stab4*stab6 - stab2*stab3*stab4*stab5),stab1*stab4 - stab2*stab3)); 
    

    std::cout << x << std::endl;

    int count = 0;
    
    for (int i = 0 ; i < 1000000 ; i ++)
    {
        Vector xrand = x.random();
        IntervalVector aux(xrand);
        if (i%10000 ==0) std::cout << i << " ; min = " << x_min << " ; max = " << x_max << std::endl;
        
        Interval res = stab->eval(aux);
        //std::cout << xrand << " " << res << std::endl;

        if (res.ub() < 0)
        {
            for (int j = 0 ; j < 5 ; j++)
            {
                if (x_min[j] > xrand[j]) x_min[j] = xrand[j];
                if (x_max[j] < xrand[j]) x_max[j] = xrand[j];
            }
            count ++;
        }
        else
        {
           //std::cout << "instable = " << xrand << std::endl; 
        }
    }

    std::cout << "count = " << count << std::endl;
}


int main(int argc, char* argv[]) {
    RNG::srand(1234);
    //test_stability();
    int count = 0;
    for (int i =1 ; i < 400000 ; i ++) {
    int p = rand();
        std::cout << p-ibex::pow(2,32) << std::endl;   
    }
    std::cout << "count = " << count << std::endl;
    return 1;
}