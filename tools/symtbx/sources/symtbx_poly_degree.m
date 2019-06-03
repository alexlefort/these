%% Give degree of a polynomial expression
%%

function d = symtbx_poly_degree(poly,p)

	[coefs,orders] = coeffs(poly,p);
	syms aux;
	orders2 = subs(poly2sym(orders,aux),{aux},{1});
	orders3 = sym2poly(orders2);
	d = length(orders3) - 1;