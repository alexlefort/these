%% Function that substitude the Laplace variable s to complex pulse iw
%% For Bode Plot

function  res = symtbx_symtf_s_to_iw(expr)

    w = sym('w','real');
    res = simplify(subs(expr, 's', 1i*w));


