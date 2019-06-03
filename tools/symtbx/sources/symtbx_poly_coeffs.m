%% Give coefficients of a polynom expression
%%

function c = symtbx_poly_coeffs(poly,p)

    [coefs,orders] = coeffs(poly,p);
	syms aux;
	orders2 = subs(poly2sym(orders,aux),{aux},{1});
	orders3 = sym2poly(orders2);
	d = length(orders3) - 1;
	orders4 =sym(orders3); 

	idx = 1;
	for ii=1:d+1
		if (orders3(ii) == 1)
			orders4(ii) = coefs(idx);
			idx = idx +1;
		end
	end
	
	c = orders4;