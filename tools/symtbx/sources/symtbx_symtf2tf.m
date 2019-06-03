%% Return a matlab transfer function from a symbolic transfert function

function tf_t = symtbx_symtf2tf(sym_tf,s)

    [num,den] = numden(sym_tf);
    
    num_c = eval(symtbx_poly_coeffs(num,s));
    den_c = eval(symtbx_poly_coeffs(den,s));
    
    if (strcmp(class(num_c),'double') ~= 1)
    	error('Remanent symbolic parameters in the symbolic transfer expression in the numerator.');
    end
    if (strcmp(class(den_c),'double') ~= 1)
    	error('Remanent symbolic parameters in the symbolic transfer expression in the denominator.');
    end
    
    tf_t = tf(num_c,den_c);