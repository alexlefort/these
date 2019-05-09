#include "ibex.h"
#include <string>
#include <ctime>
#include "ibex_Random.h"


using namespace std;
using namespace ibex;


void labrax_depth() {

    Variable x(4);
    IntervalVector x_ini(4);

    double kz0  = 0.0   ;    
    double kpi0 = 0.0   ; 
    double ipi0 = 0.0   ;   
    double kq0  = 0.0   ;
    
    double eps_ctrl = 10.0;

    x_ini[0] = Interval(kz0  - eps_ctrl , kz0  + eps_ctrl);
    x_ini[1] = Interval(kpi0 - eps_ctrl , kpi0 + eps_ctrl);
    x_ini[2] = Interval(ipi0 - eps_ctrl , ipi0 + eps_ctrl);
    x_ini[3] = Interval(kq0  - eps_ctrl , kq0  + eps_ctrl);

    double CmQ0  = -0.8116 ;
    double CzW0  = -2.6126 ;
    double CmB10 = -0.7155 ;
    double CzB10 = -1.2233 ;
    double eps   =  0.000000;

    Variable y(5);
    IntervalVector y_ini(5);
    
    y_ini[0] = Interval(CmQ0  - eps*fabs(CmQ0 ), CmQ0  + eps*fabs(CmQ0 ));
    y_ini[1] = Interval(CzW0  - eps*fabs(CzW0 ), CzW0  + eps*fabs(CzW0 ));
    y_ini[2] = Interval(CmB10 - eps*fabs(CmB10), CmB10 + eps*fabs(CmB10));
    y_ini[3] = Interval(CzB10 - eps*fabs(CzB10), CzB10 + eps*fabs(CzB10));
    y_ini[4] = Interval(-3,1);    

    Variable p(4);
    IntervalVector p_ini(4);
    
    p_ini[0] = Interval(CmQ0  - eps*fabs(CmQ0 ), CmQ0  + eps*fabs(CmQ0 ));
    p_ini[1] = Interval(CzW0  - eps*fabs(CzW0 ), CzW0  + eps*fabs(CzW0 ));
    p_ini[2] = Interval(CmB10 - eps*fabs(CmB10), CmB10 + eps*fabs(CmB10));
    p_ini[3] = Interval(CzB10 - eps*fabs(CzB10), CzB10 + eps*fabs(CzB10));

    int num_thread = 8;

    double x_prec(1e-10), y_prec(1e-10), stop_prec(0.1);

    cout << "test1" << endl;

    std::vector<Function*>         coeffs   = std::vector<Function*>();   
    std::vector<Function*>         goals    = std::vector<Function*>();   
    std::vector<Function*>         stabs    = std::vector<Function*>();   
    std::vector<NormalizedSystem*> sys_x    = std::vector<NormalizedSystem*>();
    std::vector<NormalizedSystem*> sys_xy   = std::vector<NormalizedSystem*>();
    std::vector<NormalizedSystem*> fa_y_sys = std::vector<NormalizedSystem*>();
    std::vector<Ctc*>              x_ctc_id = std::vector<Ctc*>();          
    std::vector<Ctc*>              xy_ctc   = std::vector<Ctc*>();          
    std::vector<Ctc*>              fa_y_ctc = std::vector<Ctc*>();  

    for (int i = 0 ; i <num_thread ; i++) {

        Function crit_hinf_zz1("functions/Tzz1.txt" ) ;
        Function crit_hinf_zz2("functions/Tzz2.txt" ) ;
        Function crit_hinf_zb1("functions/Tzb.txt" ) ;
    
        coeffs.push_back(new Function("functions/Tstab_coefs.txt"));

        const ExprNode& stab1 = (*coeffs[i])(x,p)[0];
        const ExprNode& stab2 = (*coeffs[i])(x,p)[1];
        const ExprNode& stab3 = (*coeffs[i])(x,p)[2];
        const ExprNode& stab4 = (*coeffs[i])(x,p)[3];
        const ExprNode& stab5 = (*coeffs[i])(x,p)[4];
        const ExprNode& stab6 = (*coeffs[i])(x,p)[5];

        stabs.push_back(new Function(x,p,ibex::max(ibex::max(ibex::max(ibex::max(-stab6,-stab4),-stab2),ibex::pow(stab1,2)*ibex::pow(stab6,2) + ibex::pow(stab2,2)*ibex::pow(stab5,2) + stab1*ibex::pow(stab4,2)*stab5 + stab2*ibex::pow(stab3,2)*stab6 - 2*stab1*stab2*stab5*stab6 - stab1*stab3*stab4*stab6 - stab2*stab3*stab4*stab5),stab1*stab4 - stab2*stab3)));
        Function f_max1(x,y,crit_hinf_zz1(x,y));
        Function f_max2(x,y,crit_hinf_zz2(x,y));
        Function f_max3(x,y,crit_hinf_zb1(x,y));
            
        goals.push_back(new Function(x,y,ibex::max(f_max1(x,y),(ibex::max(f_max2(x,y),f_max3(x,y))))));

        SystemFactory fac_x;
        fac_x.add_var(x,x_ini);
             
        SystemFactory fac_xy;
        fac_xy.add_var(x,x_ini);
        fac_xy.add_var(y,y_ini);
        fac_xy.add_goal(*goals[i]);
            
        SystemFactory fac_fa_y_sys;
        fac_fa_y_sys.add_var(x,x_ini);
        fac_fa_y_sys.add_var(p,p_ini);
        fac_fa_y_sys.add_goal(*stabs[i]);

        sys_x.push_back(     new NormalizedSystem(fac_x));       
        sys_xy.push_back(    new NormalizedSystem(fac_xy));       
        fa_y_sys.push_back(  new NormalizedSystem(fac_fa_y_sys));
    
        xy_ctc.push_back(    new CtcHC4(*sys_xy[i]));                    
        x_ctc_id.push_back(  new CtcIdentity(x_ini.size()));             
        fa_y_ctc.push_back(  new CtcIdentity(x_ini.size()+p_ini.size()));       
    }

    double prec_fa_y = 1e-5;

    OptimMinMax oo(sys_x, sys_xy, fa_y_sys, x_ctc_id, xy_ctc, fa_y_ctc, x_prec, y_prec,stop_prec,prec_fa_y, num_thread);
    
    oo.list_elem_absolute_max = 500;
    oo.list_rate = 0;
    oo.critpr = 0.4;
    oo.heap_prob = 0.5;
    oo.min_prec_coef = 10;
    oo.iter = 10;
    oo.visit_all = false;
    oo.nb_point = 1;

    oo.list_elem_absolute_max_csp = 500;
    oo.iter_csp = 10;
    oo.critpr_csp = 0.3;
    oo.list_rate_csp = 0;
    oo.min_prec_coef_csp = 10;
    oo.visit_all_csp = false;

    oo.trace=1;
    oo.trace_freq = 1;
    oo.timeout=600;

    Optim::Status res = oo.optimize(x_ini);

    oo.report();

}


int main(int argc, char* argv[]) {

RNG::srand(1234);
labrax_depth();

}
