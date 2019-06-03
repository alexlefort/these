function bool = test_poly_degree

	syms a b c d p;
	
	poly = a^3*p+2*(p^4*d-sqrt(a))*p^2*c + p^2*d^4 - 5*c*p + a*b*3;
	
	deg_res    = 6;
	
	deg_res2  = symtbx_poly_degree(poly,p);
	bool      = isequal(deg_res,deg_res2);