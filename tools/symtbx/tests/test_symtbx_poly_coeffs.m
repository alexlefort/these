function bool = test_poly_coeffs
	
	syms a b c d p;
	
	poly = a^3*p+2*(p^4*d-sqrt(a))*p^2*c + p^2*d^4 - 5*c*p + a*b*3;
	
	pol_res_2 = factor(simplify(poly));
	coefs     = symtbx_poly_coeffs(poly,p);
	pol_res   = factor(simplify(poly2sym(coefs,p)));
	bool      = isequal(pol_res,pol_res_2);