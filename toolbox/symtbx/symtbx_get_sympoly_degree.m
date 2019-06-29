function d = symtbx_get_sympoly_degree(poly,s)

    %% Give degree of a polynomial expression
    if (nargin < 2)
       syms s;
    end

    [coefs,orders] = coeffs(poly,s);
    syms aux;
    orders2 = subs(poly2sym(orders,aux),{aux},{1});
    orders3 = sym2poly(orders2);
    d = length(orders3) - 1;
