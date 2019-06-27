%% Function that substitude the complex pulse iw to Laplace variable s

function  res = symtbx_symtf_iw_to_s(expr)

    syms s;
    res = simplify(subs(expr, 'w', -1i*s));




