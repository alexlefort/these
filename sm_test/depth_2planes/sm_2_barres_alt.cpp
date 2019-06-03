#include "ibex.h"
#include <string>
#include <ctime>
#include "ibex_Random.h"


using namespace std;
using namespace ibex;

Function create_function(Function& f, Variable& x, Variable& y) {
    
    Function fres(x,y,f(x,y,Interval(-3,-3)));
    int n = 10;

    std::vector<Function> fvect = std::vector<Function>();
    fvect.push_back(Function(x,y,fres(x,y)));

    for (int i = 1 ; i <n; i ++) {
        double z = -3.0+0.4*i;
        Function new_f(x,y,f(x,y,Interval(z,z)));
        fvect.push_back(Function(x,y,ibex::max(fvect[i-1](x,y),new_f(x,y))));
    }

    return (fvect[n-1]);
}


void sm_2_barres() {

    Variable x(5);
    IntervalVector x_ini(5);
    
    double eps_ctrl = 10.0;

    x_ini[0] = Interval(-eps_ctrl,eps_ctrl);
    x_ini[1] = Interval(-eps_ctrl,eps_ctrl);
    x_ini[2] = Interval(-eps_ctrl,eps_ctrl);
    x_ini[3] = Interval(-eps_ctrl,eps_ctrl);
    x_ini[4] = Interval(-eps_ctrl,eps_ctrl);

    double CzW0      = -2.7  ;
    double CmQ0      = -0.38 ;
    double eps       =  0.0  ;
    double eps_stab  =  0.0  ;

    Variable y(2);
    IntervalVector y_ini(2);
    
    y_ini[0] = Interval(CzW0 - eps*fabs(CzW0), CzW0 + eps*fabs(CzW0));
    y_ini[1] = Interval(CmQ0 - eps*fabs(CmQ0), CmQ0 + eps*fabs(CmQ0));

    Variable p(2);
    IntervalVector p_ini(2);
    
    p_ini[0] = Interval(CzW0 - eps_stab*fabs(CzW0), CzW0 + eps_stab*fabs(CzW0));
    p_ini[1] = Interval(CmQ0 - eps_stab*fabs(CmQ0), CmQ0 + eps_stab*fabs(CmQ0));

    int num_thread = 4;

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

        Function crit_hinf_zz1("functions/Tzz1.txt" ) ;
        Function crit_hinf_zz2("functions/Tzz2.txt" ) ;
        Function crit_hinf_zb1("functions/Tzb1.txt" ) ;
        Function crit_hinf_zb2("functions/Tzb2.txt" ) ; 
    
        coeffs.push_back(new Function("functions/Tstab_coefs.txt"));

        const ExprNode& stab1 = (*coeffs[i])(x,p)[0];
        const ExprNode& stab2 = (*coeffs[i])(x,p)[1];
        const ExprNode& stab3 = (*coeffs[i])(x,p)[2];
        const ExprNode& stab4 = (*coeffs[i])(x,p)[3];
        const ExprNode& stab5 = (*coeffs[i])(x,p)[4];
        const ExprNode& stab6 = (*coeffs[i])(x,p)[5];

        stabs.push_back(new Function(x,p,ibex::max(ibex::max(ibex::max(ibex::max(-stab6,-stab4),-stab2),ibex::pow(stab1,2)*ibex::pow(stab6,2) + ibex::pow(stab2,2)*ibex::pow(stab5,2) + stab1*ibex::pow(stab4,2)*stab5 + stab2*ibex::pow(stab3,2)*stab6 - 2*stab1*stab2*stab5*stab6 - stab1*stab3*stab4*stab6 - stab2*stab3*stab4*stab5),stab1*stab4 - stab2*stab3)));

        Function fres1 = create_function(crit_hinf_zz1, x, y);
        Function fres2 = create_function(crit_hinf_zz2, x, y);
        Function fres3 = create_function(crit_hinf_zb1, x, y);
        Function fres4 = create_function(crit_hinf_zb2, x, y);

        //Function fres1(x,y,crit_hinf_zz1(x, y));
        //Function fres2(x,y,crit_hinf_zz2(x, y));
        //Function fres3(x,y,crit_hinf_zb1(x, y));
        //Function fres4(x,y,crit_hinf_zb2(x, y));

        goals.push_back(new Function(x,y,ibex::max(fres1(x,y),(ibex::max(fres2(x,y),ibex::max(fres3(x,y),fres4(x,y)))))));
        //goals.push_back(new Function(x,y,fres1(x,y)));

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
    oo.timeout=600;
    oo.eval_period  = 1;
    Optim::Status res = oo.optimize(x_ini);

    oo.report();

}


int main(int argc, char* argv[]) {

RNG::srand(1234);
sm_2_barres();

}
