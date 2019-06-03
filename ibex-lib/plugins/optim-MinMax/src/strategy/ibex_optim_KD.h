#ifndef OPTIM_KD_H
#define OPTIM_KD_H


#include "ibex_dabbene.h"
#include "ibex_kharitonov.h"
#include "ibex_Heap.h"
#include "ibex_Bsc.h"



namespace ibex {


class CostFuncKD : public CostFunc<Cell> { // element are sorted from the lowest lb of the evaluation of the objective function to the greatest
public:
        CostFuncKD();
        virtual double cost(const Cell& elem) const;
};


int KD_operator(IntervalVector& v);
int test_operator(IntervalVector& v);

class OptimKD {

private:

	int nb_iter;
	NormalizedSystem f;
    IntervalVector x_ini;
    
    Heap<Cell>* heap_stable;
    Heap<Cell>* heap_uncertain;
    Heap<Cell>* heap_unstable;
    
    Bsc* bsc;

    double volume_uncertain;
    double volume_stable;
    double volume_unstable;
    int trace_freq;

public:
	OptimKD(NormalizedSystem& f_t, IntervalVector& x_ini_t, int nb_iter_t, int trace_freq_t);
    void optimize();
    void optimize_alt();
};

}


#endif
