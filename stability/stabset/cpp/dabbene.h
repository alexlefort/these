#ifndef DABBENE_H
#define DABBENE_H


#include "ibex.h"
#include "Eigen/Eigen"
#include <linprog.h>
#include <sort.h>


bool dabbene(ibex::IntervalVector& poly, int max_iter);

ibex::IntervalVector get_odd_coeffs(ibex::IntervalVector const& poly);
ibex::IntervalVector get_even_coeffs(ibex::IntervalVector const& poly);
bool get_ko(Eigen::MatrixXd& res, ibex::IntervalVector const& odd_coeffs, int max_iter);

Eigen::MatrixXd compute_po(Eigen::MatrixXd const& ko);
Eigen::MatrixXd compute_po_roots(Eigen::MatrixXd const& p);

ibex::Interval is_interval_empty(ibex::Interval k2ip1, double k2im1, double k2im3, double i, double no);

Eigen::MatrixXd compute_ve(Eigen::MatrixXd const & odd_roots, int ne);
Eigen::MatrixXd compute_u(double omega, int n);

bool does_ke_exist(ibex::IntervalVector const& even_coeffs, Eigen::MatrixXd const& ve);


Eigen::VectorXcd compute_poly_roots(Eigen::MatrixXd const& p);

// Test if values in a vector a real
bool test_if_real(Eigen::VectorXcd const& p);
// Sort values in a vector
Eigen::MatrixXd sort_vector(Eigen::MatrixXd const& p);
bool test_if_distinct(Eigen::MatrixXd const& p);


Eigen::MatrixXd return_roots(Eigen::MatrixXd const& r);


#endif