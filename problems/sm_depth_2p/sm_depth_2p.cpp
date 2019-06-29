#include "ibex.h"
#include <string>
#include <ctime>
#include "ibex_Random.h"


using namespace std;
using namespace ibex;

void sm_2_barres() {

    Variable x(9);
    IntervalVector x_ini(9);

    double eps_ctrl = 30;

    x_ini[0] = Interval(-eps_ctrl,eps_ctrl);
    x_ini[1] = Interval(-eps_ctrl,eps_ctrl);
    x_ini[2] = Interval(-eps_ctrl,eps_ctrl);
    x_ini[3] = Interval(-eps_ctrl,eps_ctrl);
    x_ini[4] = Interval(-eps_ctrl,eps_ctrl);
    x_ini[5] = Interval(-eps_ctrl,eps_ctrl);
    x_ini[6] = Interval(-eps_ctrl,eps_ctrl);
    x_ini[7] = Interval(-eps_ctrl,eps_ctrl);
    x_ini[8] = Interval(-0.01,0.01);
    
    double CmQ0       = -0.38   ;
    double CzW0       = -2.7    ;
    double CmB10      = -0.3230 ;
    double CmB20      =  0.0960 ;
    double CzB10      = -0.7030 ;
    double CzB20      = -0.5040 ;

    double eps        =  0.0  ;
    double eps2       =  0.02 ;
    Variable y(7);
    IntervalVector y_ini(7);
    
    y_ini[0] = Interval(CmQ0  - eps*fabs(CmQ0 ), CmQ0  + eps*fabs(CmQ0 ));
    y_ini[1] = Interval(CzW0  - eps*fabs(CzW0 ), CzW0  + eps*fabs(CzW0 ));
    y_ini[2] = Interval(CmB10 - eps*fabs(CmB10), CmB10 + eps*fabs(CmB10));
    y_ini[3] = Interval(CmB20 - eps*fabs(CmB20), CmB20 + eps*fabs(CmB20));
    y_ini[4] = Interval(CzB10 - eps*fabs(CzB10), CzB10 + eps*fabs(CzB10));
    y_ini[5] = Interval(CzB20 - eps*fabs(CzB20), CzB20 + eps*fabs(CzB20));
    y_ini[6] = Interval(-3,1);

    Variable p(6);
    IntervalVector p_ini(6);
    
    p_ini[0] = Interval(CmQ0  - eps2*fabs(CmQ0 ), CmQ0  + eps2*fabs(CmQ0 ));
    p_ini[1] = Interval(CzW0  - eps2*fabs(CzW0 ), CzW0  + eps2*fabs(CzW0 ));
    p_ini[2] = Interval(CmB10 - eps2*fabs(CmB10), CmB10 + eps2*fabs(CmB10));
    p_ini[3] = Interval(CmB20 - eps2*fabs(CmB20), CmB20 + eps2*fabs(CmB20));
    p_ini[4] = Interval(CzB10 - eps2*fabs(CzB10), CzB10 + eps2*fabs(CzB10));
    p_ini[5] = Interval(CzB20 - eps2*fabs(CzB20), CzB20 + eps2*fabs(CzB20));

    int num_thread = 8;

    double x_prec(1e-6), y_prec(1e-6), stop_prec(0.1);

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

    cout << "test2" << endl;      

    for (int i = 0 ; i <num_thread ; i++) {

        Function zzc("functions/Tzzc.txt" ) ;
        Function zzs("functions/Tzzs.txt" ) ; 
        Function ttc("functions/Tttc.txt" ) ; 

        Function tb1("functions/Ttb1c.txt" ) ;
        Function tb2("functions/Ttb2c.txt" ) ;
        Function zb1("functions/Tzb1c.txt" ) ;
        Function zb2("functions/Tzb2c.txt" ) ;

        Function z(x,y,ibex::max(zzc(x,y),(ibex::max(zzs(x,y),ibex::max(zb1(x,y),zb2(x,y))))));
        Function t(x,y,ibex::max(ttc(x,y),(ibex::max(tb1(x,y),tb2(x,y)))));

        coeffs.push_back(new Function("functions/Tstab_coefs.txt"));
            
        goals.push_back(new Function(x,y,ibex::max(z(x,y),t(x,y))));

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

    cout << "start" << endl;

    OptimMinMax oo(sys_x, sys_xy, fa_y_sys, x_ctc_id, xy_ctc, fa_y_ctc, x_prec, y_prec,stop_prec,prec_fa_y, num_thread);
    
    oo.list_elem_absolute_max = 0;
    oo.list_rate = 0;
    oo.critpr = 0.0;
    oo.heap_prob = 0;
    oo.min_prec_coef = 1;
    oo.iter = 10;
    oo.visit_all = true;
    oo.nb_point = 1;

    oo.list_elem_absolute_max_csp = 0;
    oo.iter_csp = 10;
    oo.critpr_csp = 0.0;
    oo.list_rate_csp = 0;
    oo.min_prec_coef_csp = 1;
    oo.visit_all_csp = false;

    oo.trace=1;
    oo.trace_freq = 10;
    oo.timeout=1500;
    oo.eval_period  = 1;
    Optim::Status res = oo.optimize(x_ini);

    oo.report();

}


int main(int argc, char* argv[]) {

RNG::srand(1234);
sm_2_barres();

}
