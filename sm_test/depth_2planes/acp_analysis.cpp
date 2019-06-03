#include "ibex.h"
#include <string>
#include <ctime>
#include "ibex_Random.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace ibex;


void acp_analysis() {

    ofstream myfile;
    myfile.open ("acp_result.txt");
    

    Variable p(11);
    IntervalVector x(11);

    double xminval = -10.0;
    double xmaxval =  10.0;

    for (int i = 0 ; i < 9 ; i++)
    {
    	x[i] = Interval(xminval, 3*xmaxval);
    }

    //double kpz10 = 1.0712    ;    
    //double kdz10 = 0.9998    ; 
    //double kpt10 = -10.3220  ;   
    //double kdt10 = -9.0197   ;
    //double kiz10 =   0.001   ;  
    //double kpz20 = -0.3585   ;  
    //double kdz20 = -0.0221   ;  
    //double kpt20 = -1.8123   ;  
    //double kdt20 = -0.8398   ;  
    //
    double eps_ctrl = 0.0;
//
    //x[0] = Interval(kpz10 - eps_ctrl*fabs(kpz10), kpz10 + 3*eps_ctrl*fabs(kpz10));
    //x[1] = Interval(kdz10 - eps_ctrl*fabs(kdz10), kdz10 + 3*eps_ctrl*fabs(kdz10));
    //x[2] = Interval(kpt10 - eps_ctrl*fabs(kpt10), kpt10 + 3*eps_ctrl*fabs(kpt10));
    //x[3] = Interval(kdt10 - eps_ctrl*fabs(kdt10), kdt10 + 3*eps_ctrl*fabs(kdt10));
    //x[4] = Interval(kiz10 - eps_ctrl*fabs(kiz10), kiz10 + 3*eps_ctrl*fabs(kiz10));
    //x[5] = Interval(kpz20 - eps_ctrl*fabs(kpz20), kpz20 + 3*eps_ctrl*fabs(kpz20));
    //x[6] = Interval(kdz20 - eps_ctrl*fabs(kdz20), kdz20 + 3*eps_ctrl*fabs(kdz20));
    //x[7] = Interval(kpt20 - eps_ctrl*fabs(kpt20), kpt20 + 3*eps_ctrl*fabs(kpt20));
    //x[8] = Interval(kdt20 - eps_ctrl*fabs(kdt20), kdt20 + 3*eps_ctrl*fabs(kdt20));

    double CzW0 = -2.7  ;
    double CmQ0 = -0.38 ;
    
    x[9]  = Interval(CzW0, CzW0);
    x[10] = Interval(CmQ0, CmQ0);
   
    Function coeffs =  Function("functions/Tstab_coefs_x.txt");

    const ExprNode& stab1 = coeffs(p)[0];
    const ExprNode& stab2 = coeffs(p)[1];
    const ExprNode& stab3 = coeffs(p)[2];
    const ExprNode& stab4 = coeffs(p)[3];
    const ExprNode& stab5 = coeffs(p)[4];
    const ExprNode& stab6 = coeffs(p)[5];

    Vector x_min(9), x_max(9);

    for (int i = 0 ; i < 9 ; i++)
    {
    	x_min[i] = xmaxval;
    	x_max[i] = xminval;
    }


    Function stab = Function(p,ibex::max(ibex::max(ibex::max(ibex::max(-stab6,-stab4),-stab2),ibex::pow(stab1,2)*ibex::pow(stab6,2) + ibex::pow(stab2,2)*ibex::pow(stab5,2) + stab1*ibex::pow(stab4,2)*stab5 + stab2*ibex::pow(stab3,2)*stab6 - 2*stab1*stab2*stab5*stab6 - stab1*stab3*stab4*stab6 - stab2*stab3*stab4*stab5),stab1*stab4 - stab2*stab3));

    int count = 0;
    
    for (int i = 0 ; i < 1000000 ; i ++)
    {
        Vector xrand = x.random();
        IntervalVector aux(xrand);
        if (i%10000 ==0) std::cout << i << " ; min = " << x_min << " ; max = " << x_max << std::endl;
        
        Interval res = stab.eval(aux);
        //std::cout << xrand << " " << res << std::endl;

        if (res.ub() < 0)
        {
            for (int j = 0 ; j < 9 ; j++)
            {
                if (x_min[j] > xrand[j]) x_min[j] = xrand[j];
                if (x_max[j] < xrand[j]) x_max[j] = xrand[j];
                myfile << xrand[j] << ",";
            }
            count ++;

            myfile << " \n" ;
    
        }
        else
        {
           //std::cout << "instable = " << xrand << std::endl; 
        }
    }

    std::cout << "count = " << count << std::endl;
    myfile.close();
}


int main(int argc, char* argv[]) {

RNG::srand(1234);
acp_analysis();

}
