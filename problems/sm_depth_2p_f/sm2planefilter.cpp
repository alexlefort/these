#include "ibex.h"
#include <string>
#include <ctime>
#include "ibex_Random.h"


using namespace std;
using namespace ibex;


void sm_2_barres() {

    Variable x(13);
    IntervalVector x_ini(13);
    
    double eps_ctrl = 10.0;

    for (int i = 0 ; i < 13 ; i ++)
    {
        x_ini[i] = Interval(-eps_ctrl, eps_ctrl);
    }
    x_ini[5]  = Interval(-1,1);
    x_ini[9]  = Interval(1e-5,30);
    x_ini[10] = Interval(1e-5,30);
    x_ini[11] = Interval(0.9,1.0);
    x_ini[12] = Interval(0.9,1.0);

    double CzW0      = -2.7  ;
    double CmQ0      = -0.38 ;
    double eps       =  0.0  ;
    double eps_stab  =  0.0  ;

    Variable y(3);
    IntervalVector y_ini(3);
    
    y_ini[0] = Interval(CzW0 - eps*fabs(CzW0), CzW0 + eps*fabs(CzW0));
    y_ini[1] = Interval(CmQ0 - eps*fabs(CmQ0), CmQ0 + eps*fabs(CmQ0));
    y_ini[2] = Interval(-2,1);

    Variable p(2);
    IntervalVector p_ini(2);
    
    p_ini[0] = Interval(CzW0 - eps_stab*fabs(CzW0), CzW0 + eps_stab*fabs(CzW0));
    p_ini[1] = Interval(CmQ0 - eps_stab*fabs(CmQ0), CmQ0 + eps_stab*fabs(CmQ0));

    int num_thread = 4;

    double x_prec(1e-6), y_prec(1e-6), stop_prec(0.1);

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

        Function zz1("functions/Tzz1.txt" ) ;
        Function zz2("functions/Tzz2.txt" ) ;
        Function zb1("functions/Tzb1.txt" ) ;
        Function zb2("functions/Tzb2.txt" ) ; 
    
        coeffs.push_back(new Function("functions/Tstab_coefs.txt"));


        goals.push_back(new Function(x,y,ibex::max(zz1(x,y),(ibex::max(zz2(x,y),ibex::max(zb1(x,y),zb2(x,y)))))));

        SystemFactory fac_x;
        fac_x.add_var(x,x_ini);
             
        SystemFactory fac_xy;
        fac_xy.add_var(x,x_ini);
        fac_xy.add_var(y,y_ini);
        fac_xy.add_goal(*goals[i]);

        SystemFactory fac_fa_y_sys;
        fac_fa_y_sys.add_var(x,x_ini);
        fac_fa_y_sys.add_var(p,p_ini);
        fac_fa_y_sys.add_goal(*coeffs[i]);

        sys_x.push_back(     new NormalizedSystem(fac_x));       
        sys_xy.push_back(    new NormalizedSystem(fac_xy));       
        fa_y_sys.push_back(  new NormalizedSystem(fac_fa_y_sys));
    
        xy_ctc.push_back(    new CtcHC4(*sys_xy[i]));                    
        x_ctc_id.push_back(  new CtcIdentity(x_ini.size()));             
        fa_y_ctc.push_back(  new CtcIdentity(x_ini.size()+p_ini.size()));       
    }

    double prec_fa_y = 1e-8;

    OptimMinMax oo(sys_x, sys_xy, fa_y_sys, x_ctc_id, xy_ctc, fa_y_ctc, x_prec, y_prec,stop_prec,prec_fa_y, num_thread);
    
    oo.list_elem_absolute_max = 0;
    oo.list_rate = 0;
    oo.critpr = 0.1;
    oo.heap_prob = 0;
    oo.min_prec_coef = 10;
    oo.iter = 10;
    oo.visit_all = true;
    oo.nb_point = 1;

    oo.list_elem_absolute_max_csp = 0;
    oo.iter_csp = 10;
    oo.critpr_csp = 0.1;
    oo.list_rate_csp = 0;
    oo.min_prec_coef_csp = 10;
    oo.visit_all_csp = true;

    oo.trace=1;
    oo.trace_freq = 10;
    oo.timeout=1000;
    oo.eval_period  = 1;
    Optim::Status res = oo.optimize(x_ini);

    oo.report();

}


int main(int argc, char* argv[]) {

RNG::srand(1234);
sm_2_barres();

}
