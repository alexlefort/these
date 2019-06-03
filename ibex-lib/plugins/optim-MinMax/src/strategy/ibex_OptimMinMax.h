//============================================================================
//                                  I B E X
// File        : ibex_OptimMinMax.h
// Author      : Dominique Monnet, Jordan Ninin
// License     : See the LICENSE file
// Created     : Oct 1, 2016
//============================================================================


#ifndef __OPTIM_MINIMAX__
#define __OPTIM_MINIMAX__


#include "ibex_Ctc.h"
#include "ibex_IntervalVector.h"
#include "ibex_NormalizedSystem.h"
#include "ibex_LightOptimMinMax.h"
#include "ibex_Bsc.h"
#include "ibex_Optim.h"
#include "omp.h"

namespace ibex {

class OptimMinMax : public Optim {

public:

    /* Constructor*/
    OptimMinMax(std::vector<NormalizedSystem*>  x_sys_t        , 
                std::vector<NormalizedSystem*>  xy_sys_t       , 
                std::vector<Ctc*>               x_ctc_t        , 
                std::vector<Ctc*>               xy_ctc_t       ,
                double                          prec_x         , 
                double                          prec_y         , 
                double                          goal_rel_prec  ,
                int                             num_thread_t   );

    /* Constructor, with "for all y" constraints
     * for all constraint in the objectif function of max_fa_y_cst lower than 0. If several for all constraints the objectif is equal to the max
     * of constraint functions lower than 0*/
    OptimMinMax(std::vector<NormalizedSystem*> x_sys_t            , 
                std::vector<NormalizedSystem*> xy_sys_t           , 
                std::vector<NormalizedSystem*> max_fa_y_cst_sys_t , 
                std::vector<Ctc*>              x_ctc_t            , 
                std::vector<Ctc*>              xy_ctc_t           , 
                std::vector<Ctc*>              y_fa_ctc_t         ,
                double                         prec_x             , 
                double                         prec_y             , 
                double                         goal_rel_prec      , 
                double                         fa_cst_prec        ,
                int                            num_thread_t      );

    /* Runs a B&B like algorithm
     * arguments: -x_ini: initial x box
     *            -y_ini: initial y box
     *            -prec_x: minimum size of x box
     *            -prec_y: minimum size of y box
     *            -stop_prec: stop criterion
     * */

	/**
	 * \brief Run the optimization.
	 *
	 * \param init_box             The initial box
	 * \param obj_init_bound       (optional) can be set when an initial upper bound of the objective minimum is known a priori.
	 *                             This bound can be obtained, e.g., by a local solver. This is equivalent to (but more practical
	 *                             than) adding a constraint f(x)<=obj_init_bound.
	 *
	 * \return SUCCESS             If the global minimum (with respect to the precision required) has been found.
	 *                             In particular, at least one feasible point has been found, less than obj_init_bound, and in the time limit.
	 *
	 *         INFEASIBLE          if no feasible point exist less than obj_init_bound. In particular, the function returns INFEASIBLE
	 *                             if the initial bound "obj_init_bound" is LESS than the true minimum (this case is only possible if
	 *                             goal_abs_prec and goal_rel_prec are 0). In the latter case, there may exist feasible points.
	 *
	 *         NO_FEASIBLE_FOUND   if no feasible point could be found less than obj_init_bound. Contrary to INFEASIBLE,
	 *                             infeasibility is not proven here. Warning: this return value is sensitive to the goal_abs_prec and
	 *                             goal_rel_prec parameters. The upperbounding makes the optimizer only looking for points less than
	 *                             min { (1-goal_rel_prec)*obj_init_bound, obj_init_bound - goal_abs_prec }.
	 *
	 *         UNBOUNDED_OBJ       the objective function seems unbounded (tends to -oo).
	 *
	 *         TIMEOUT             time is out.
	 */

    virtual ~OptimMinMax(); // Destructor

    IntervalVector x_box_init         ; // x box where solution is searched
    IntervalVector y_box_init         ; // y box condition
    IntervalVector y_box_init_fa      ; // y box for constraints of type for all y may be different than y_box_ini

    bool   propag                     ; // true if propagation use
    int    trace_freq                 ;
    int    num_thread                 ;

                                        //light solver param for objectif function
    int    list_rate                  ; // rate of growth of y_heap size, see compute_heap_max_size function for formula detail
    int    list_elem_absolute_max     ;
    int    iter                       ; // number of iteration of lightsolver allowedNormalizedSystem *
    double min_prec_coef              ; // used to compute y_prec allowed when run the light solver, see compute_min_prec function for formula
    int    critpr                     ; // probability to choose second heap in light_solver

                                        //light solver param for csp
    int    list_rate_csp              ; // rate of growth of y_heap size, see compute_heap_max_size function for formula detail
    int    list_elem_absolute_max_csp ;
    int    iter_csp                   ; // number of iteration of lightsolver allowedNormalizedSystem *
    double min_prec_coef_csp          ; // used to compute y_prec allowed when run the light solver, see compute_min_prec function for formula
    int    critpr_csp                 ; // probability to choose second heap in light_solver
    bool   only_csp                   ; // solve a csp problem, no objective function.

    int    heap_prob                  ; // probability between 0 and 100 to choose heap 2 over 1 in the main B&B, heap 1 is sorted from the lower lb
                                        // to the greatest, heap 2 from the greatest ub to the lowest
    bool   visit_all                  ; // if true force to visit all leaves of slave heap in max optim
    bool   visit_all_csp              ; // if true force to visit all leaves of slave heap in max optim of csp solver

    int    nb_point                   ; // number of point evaluation in current box x (use to compute uplo)
    int eval_period                   ; // period of point evaluation performed (expected to be greater ore equal to one)
    int eval_counter                  ; // counter linked to evaluation period
    double perf_thresh                ; // cell is spawn from loc_solve csp list if its perf_csp is lower than perf_thresh

    double prec_y                     ;
    bool   fa_y_cst                   ; // indicates if there is for all y constraints
    bool   min_goal                   ; // true if minimization problem, false if min max problem
    double prec_fa_y                  ;

    omp_lock_t bufferlock ;

//private:
    std::vector<Ctc*>               x_ctc     ; // contractor w.r.t constraint on x
    std::vector<NormalizedSystem*>  x_sys     ; // contains cst on x and objective function
    std::vector<LightOptimMinMax*>  lsolve    ;

    std::vector<Bsc*>               bsc       ;
    std::vector<Function*>          minus_goal_y_at_x ; // goal function f becomes -f to solve a minimization problem over y at a fixed x

    // Cst Factory
    std::vector<LightOptimMinMax*>  fa_lsolve              ;
    std::vector<Function*>          minus_goal_csp_y_at_x  ; // goal csp function g becomes -g to solve a minimization problem over y at a fixed x

    Status optimize(const IntervalVector& init_box, double obj_init_bound=POS_INFINITY);
    Status optimize();

    std::vector<Cell*> nsect_cell(int n, Cell *c);
    inline void set_prec_y(double prec_y) {this->prec_y = prec_y; }

    void   init_lsolve();
    void   init_fa_lsolve();
    double compute_min_prec( const IntervalVector& x_box,bool csp);
    int    choose_nbiter(bool midpoint_eval,bool csp,Cell* x_cell);
    int    compute_heap_max_size(int y_heap_size,bool csp);
    bool   get_feasible_point(Cell * elem, int ith);
    int    check_constraints(Cell * x_cell,bool midp, int ith);
    int    check_regular_ctr(const IntervalVector& box, int ith);
    int    check_fa_ctr(Cell * x_cell,bool midp, int ith);
    bool   handle_cell(Cell * x_cell, int ith);
    bool   spawn(Cell* x_cell); // spawn a cell into subcells from the csp list of loc_solve if perf crit is lower than spawn_thresh
    bool   check_optimizer();
    bool   eval_goal(Cell* x_cell,double loup, int ith);
    void   show_yheap(Cell * x_cell, int ith);

    //Default parameters for light optim min max solver
    static const int    default_list_rate;
    static const int    default_list_elem_absolute_max;
    static const int    default_iter;
    static const double default_min_prec_coef;
    static const int    default_prob_heap;
    static const bool   default_visit_all;
    static const int    default_nb_point;
    static const int    default_eval_period;
    static const double default_perf_thresh;

    //Csp fa function default parameters for light solver
    static const int    default_list_rate_csp;
    static const int    default_list_elem_absolute_max_csp;
    static const int    default_iter_csp;
    static const double default_min_prec_coef_csp;
    static const int    default_prob_heap_csp;
    static const bool   default_visit_all_csp;
   double volume;
};


}// end namespace

#endif
