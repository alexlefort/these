#include "optim_KD.h"


using namespace ibex;
using namespace std;

CostFuncKD::CostFuncKD() {}

double CostFuncKD::cost(const Cell& elem) const {
    return -elem.box.volume();
}

int ibex::KD_operator(IntervalVector& v) {

    bool res_k = kharitonov(v);
    bool res_d = dabbene(v,20);
    if (res_k) return 2; // Completely stable
    if (res_d) {
    	return 1; // Uncertain
    } else {
    	return 0; // Completely unstable
    }
}

int ibex::test_operator(IntervalVector& v) {
    if (v[0].lb() > 0) {
    	return 0;
    } else if (v[0].ub() < 0) {
    	return 2;
    } else {
    	return 1;
    }
}

OptimKD::OptimKD(NormalizedSystem& f_t,
	             IntervalVector&   x_ini_t,
	             int               nb_iter_t,
	             int               trace_freq_t) :
f(f_t),
x_ini(x_ini_t),
nb_iter(nb_iter_t),
trace_freq(trace_freq_t),
heap_stable(   new Heap<Cell>(* new CostFuncKD())),
heap_uncertain(new Heap<Cell>(* new CostFuncKD())),
heap_unstable( new Heap<Cell>(* new CostFuncKD()))
{
    this->bsc = new LargestFirst()  ;

    volume_uncertain = 0 ;
    volume_stable    = 0 ;
    volume_unstable  = 0 ;
    cout << "build" << endl;
}

void OptimKD::optimize_alt() {

// Initialization
	heap_uncertain->flush();
	heap_stable->flush();
	heap_unstable->flush();

	Cell x_ini_cell(x_ini);
	heap_uncertain->push(&x_ini_cell);
	volume_uncertain += x_ini.volume();

	cout << x_ini << endl;
	cout << x_ini.volume() << endl;

    //cout << "uncertain = " << volume_uncertain <<  " ; stable = " << volume_stable << " ; instable = " << volume_unstable << endl;
        
    //cout << "begin loop" << endl;
	for (int i = 0 ; i < nb_iter; i ++)
	{
		Cell* v = heap_uncertain->pop();
		volume_uncertain -= v->box.volume();
        
        //cout << "s1" << endl;

		pair<Cell*,Cell*> x = bsc->bisect_cell(*v);
		
		//cout << "s2" << endl;

		//cout << "s3" << endl;

		 IntervalVector x1 = x.first->box;
		 IntervalVector x2 = x.second->box;
// 
		 //cout << "s4" << endl;
		 IntervalVector res1 = f.goal->eval_vector(x.first->box);
		 IntervalVector res2 = f.goal->eval_vector(x.second->box);
        
        //cout << "s5 " << res1 << " " << res2 << endl;

        int res_kd_1 = test_operator(res1);
        int res_kd_2 = test_operator(res2);
        //cout << "s6 "  << res_kd_1 << " " << res_kd_2 << endl;

        if (res_kd_1 == 0) {
        	heap_unstable->push(x.first);
        	volume_unstable += x1.volume();
        } 
        if (res_kd_1 == 1) {
        	heap_uncertain->push(x.first);
        	volume_uncertain += x1.volume();
        }
        if (res_kd_1 == 2) {
        	heap_stable->push(x.first);
        	volume_stable += x1.volume();
        }
// 
        if (res_kd_2 == 0) {
        	heap_unstable->push(x.second);
        	volume_unstable += x2.volume();
        } 
        if (res_kd_2 == 1) {
        	heap_uncertain->push(x.second);
        	volume_uncertain += x2.volume();
        }
        if (res_kd_2 == 2) {
        	heap_stable->push(x.second);
        	volume_stable += x2.volume();
        }
//      
        //cout << "s7" << endl;

        if (i%trace_freq == 0){
        	cout << "uncertain = " << volume_uncertain <<  " ; stable = " << volume_stable << " ; instable = " << volume_unstable << endl;
        }

	}
}

void OptimKD::optimize() {

// Initialization
	heap_uncertain->flush();
	heap_stable->flush();
	heap_unstable->flush();

	Cell x_ini_cell(x_ini);
	heap_uncertain->push(&x_ini_cell);
	volume_uncertain += x_ini.volume();

	cout << x_ini << endl;
	cout << x_ini.volume() << endl;

    //cout << "uncertain = " << volume_uncertain <<  " ; stable = " << volume_stable << " ; instable = " << volume_unstable << endl;
        
    //cout << "begin loop" << endl;
	for (int i = 0 ; i < nb_iter; i ++)
	{
		Cell* v = heap_uncertain->pop();
		volume_uncertain -= v->box.volume();
        
        //cout << "s1" << endl;

		pair<Cell*,Cell*> x = bsc->bisect_cell(*v);
		
		//cout << "s2" << endl;

		//cout << "s3" << endl;

		 IntervalVector x1 = x.first->box;
		 IntervalVector x2 = x.second->box;
// 
		 //cout << "s4" << endl;
		 IntervalVector res1 = f.goal->eval_vector(x.first->box);
		 IntervalVector res2 = f.goal->eval_vector(x.second->box);
        
        //cout << "s5 " << res1 << " " << res2 << endl;

        int res_kd_1 = KD_operator(res1);
        int res_kd_2 = KD_operator(res2);
        //cout << "s6 "  << res_kd_1 << " " << res_kd_2 << endl;

        if (res_kd_1 == 0) {
        	heap_unstable->push(x.first);
        	volume_unstable += x1.volume();
        } 
        if (res_kd_1 == 1) {
        	heap_uncertain->push(x.first);
        	volume_uncertain += x1.volume();
        }
        if (res_kd_1 == 2) {
        	heap_stable->push(x.first);
        	volume_stable += x1.volume();
        }
// 
        if (res_kd_2 == 0) {
        	heap_unstable->push(x.second);
        	volume_unstable += x2.volume();
        } 
        if (res_kd_2 == 1) {
        	heap_uncertain->push(x.second);
        	volume_uncertain += x2.volume();
        }
        if (res_kd_2 == 2) {
        	heap_stable->push(x.second);
        	volume_stable += x2.volume();
        }
//      
        //cout << "s7" << endl;

        if (i%trace_freq == 0){
        	cout << "uncertain = " << volume_uncertain <<  " ; stable = " << volume_stable << " ; instable = " << volume_unstable << endl;
        }

	}
}